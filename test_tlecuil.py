"""
Verification and benchmark tests for TLECUIL.

    python -m pytest -q          (or:  python test_tlecuil.py)

Verification asks "is the math coded correctly?" and checks the code
against things known independently of it: hand calculations, conservation
laws, closed-form limits, and numerical integration.
Benchmark tests ask "does the model reproduce published designs within its
stated accuracy?" (see benchmarks.py and Theory Sec. 11 [sec:bench]).
"""
import math
from dataclasses import replace

from scipy.integrate import quad

from tlecuil import (FSPInputs, size_fsp, radiating_area_m2, optimal_cold_temperature,
                     analytic_optimum_ratio, STEFAN_BOLTZMANN)
import benchmarks

IDEAL = dict(panel_effectiveness=1.0, area_margin=0.0, hx_delta_t_k=0.0, radiator_delta_t_k=0.0)


# ---------------------------- verification ---------------------------------

def test_hand_calculation():
    # 75 kWe, eta 27 %, PMAD 90 %, isothermal 450 K, sink 210 K, eps 0.85, ideal panel.
    inp = FSPInputs(delivered_power_kwe=75, efficiency_mode="fixed", conversion_efficiency=0.27,
                    cold_temperature_k=450, sink_temperature_k=210, **IDEAL)
    r = size_fsp(inp)
    assert math.isclose(r.reactor_thermal_kwt, 75 / 0.9 / 0.27, rel_tol=1e-12)
    flux = 0.85 * STEFAN_BOLTZMANN * (450 ** 4 - 210 ** 4)          # W/m^2 per face
    assert math.isclose(r.radiator_area_m2, r.converter_heat_rejection_kwt * 1e3 / flux, rel_tol=1e-12)
    assert abs(r.reactor_thermal_kwt - 308.6) < 0.05
    assert abs(r.radiator_area_m2 - 119.7) < 0.05                    # radiating area, all faces


def test_two_sided_panel_halves_planform_not_radiating_area():
    # Regression test for the v0.1 error: the number of radiating sides must not
    # change the total radiating area; it only sets the panel planform.
    one = size_fsp(FSPInputs(radiating_sides=1))
    two = size_fsp(FSPInputs(radiating_sides=2))
    assert math.isclose(one.radiator_area_m2, two.radiator_area_m2, rel_tol=1e-12)
    assert math.isclose(two.radiator_planform_m2, two.radiator_area_m2 / 2, rel_tol=1e-12)


def test_energy_conserved_everywhere():
    for p in (1, 40, 100, 5000):
        r = size_fsp(FSPInputs(delivered_power_kwe=p, parasitic_fraction=0.07,
                               fixed_parasitic_kwe=1.5, reactor_heat_loss_fraction=0.12))
        assert abs(r.energy_balance_error_kw) < 1e-9 * p


def test_linear_scaling_without_fixed_loads():
    a40 = size_fsp(FSPInputs(delivered_power_kwe=40)).radiator_area_m2
    a400 = size_fsp(FSPInputs(delivered_power_kwe=400)).radiator_area_m2
    assert math.isclose(a400, 10 * a40, rel_tol=1e-12)


def test_fixed_load_breaks_linearity_in_the_expected_direction():
    base = FSPInputs(fixed_parasitic_kwe=2.0)
    r10 = size_fsp(replace(base, delivered_power_kwe=10))
    r100 = size_fsp(replace(base, delivered_power_kwe=100))
    assert r10.reactor_thermal_kwt * 10 > r100.reactor_thermal_kwt   # small systems pay more per kWe


def test_fluid_loop_matches_quadrature():
    q, t_in, t_out, ts, eps, eta_p = 2.0e5, 440.0, 360.0, 210.0, 0.85, 0.8
    k = eta_p * eps * STEFAN_BOLTZMANN
    numeric = (q / (t_in - t_out)) * quad(lambda t: 1 / (t ** 4 - ts ** 4), t_out, t_in)[0] / k
    assert math.isclose(radiating_area_m2(q, t_in, t_out, ts, eps, eta_p), numeric, rel_tol=1e-9)


def test_fluid_loop_limit_is_isothermal():
    iso = radiating_area_m2(1e5, 400, 400, 200, 0.85, 1)
    near = radiating_area_m2(1e5, 400, 399.99, 200, 0.85, 1)
    assert math.isclose(iso, near, rel_tol=1e-4)


def test_zero_kelvin_sink_branch():
    a = radiating_area_m2(1e5, 400, 300, 0.0, 1, 1)
    numeric = (1e5 / 100) * quad(lambda t: t ** -4, 300, 400)[0] / STEFAN_BOLTZMANN
    assert math.isclose(a, numeric, rel_tol=1e-9)


def test_small_and_large_sinks_match_quadrature():
    # Review finding: the closed form loses precision for small sinks. The
    # series branch must agree with quadrature from ~0 K up to near T_out.
    t_in, t_out = 400.0, 300.0
    for ts in (1e-9, 1e-6, 1e-3, 1e-2, 1.0, 100.0, 149.0, 151.0, 250.0, 290.0):
        numeric = (1e5 / (t_in - t_out)) * quad(lambda t: 1 / (t ** 4 - ts ** 4), t_out, t_in,
                                                epsabs=0, epsrel=1e-12)[0] / STEFAN_BOLTZMANN
        assert math.isclose(radiating_area_m2(1e5, t_in, t_out, ts, 1, 1), numeric, rel_tol=1e-9), ts


def test_optimum_matches_closed_form():
    for f in (0.3, 0.5, 0.7, 1.0):
        inp = FSPInputs(fraction_of_carnot=f, sink_temperature_k=0.0, hot_temperature_k=1000.0, **IDEAL)
        tc, _ = optimal_cold_temperature(inp)
        assert abs(tc / 1000.0 - analytic_optimum_ratio(f)) < 1e-5
    assert math.isclose(analytic_optimum_ratio(1.0), 0.75)


def test_second_law_guard():
    try:
        size_fsp(FSPInputs(efficiency_mode="fixed", conversion_efficiency=0.6,
                           hot_temperature_k=900, cold_temperature_k=450))
    except ValueError:
        return
    raise AssertionError("an efficiency above Carnot was accepted")


def test_fuel_rule_of_thumb():
    # ~1.05 g of U-235 fissioned per MWd at 200 MeV per fission
    r = size_fsp(FSPInputs(full_power_years=1))
    mwd = r.reactor_thermal_kwt / 1e3 * 365.25
    assert abs(r.u235_fissioned_kg * 1e3 / mwd - 1.05) < 0.02


# ------------------------------ benchmarks ---------------------------------

RES = benchmarks.run()


def test_calibration_reproduces_b1_exactly():
    b1 = RES["radiator"][0]
    assert b1["role"] == "calibration" and abs(b1["err_pred"]) < 1e-9


def test_independent_nasa_2025_designs_within_tolerance():
    for r in RES["radiator"]:
        if r["key"] in ("B2", "B3"):
            assert abs(r["err_pred"]) <= benchmarks.TOL_RADIATOR, r


def test_default_effectiveness_is_the_b1_calibration():
    assert abs(FSPInputs().panel_effectiveness - RES["eta_p_calibrated"]) < 0.005


def test_electrical_chain_and_first_law_within_5_percent():
    for r in RES["system"]:
        assert r["pass"], r


if __name__ == "__main__":
    import sys
    fails = 0
    for name, fn in [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]:
        try:
            fn()
            print("PASS", name)
        except AssertionError as e:
            fails += 1
            print("FAIL", name, e)
    sys.exit(1 if fails else 0)
