"""
TLECUIL -- A First-Order Fission Surface Power Sizer
=====================================================

*Tlecuilli* (Nahuatl): the three-stone hearth -- "fire pit with three stones
to support cooking vessels, hearth" (Karttunen 1992, p. 307). The three
stones stand for the three subsystems this code balances: the reactor
(fire), power conversion (the vessel it heats), and heat rejection (the
heat that must leave).

TLECUIL is an open, independent, first-order steady-state thermal--electrical
balance for a fission surface power (FSP) system. Every equation is derived
in docs/tlecuil_theory.pdf; comments cite the theory equation number and the
primary source as [BibTeX key], which match docs/references.bib.

It is NOT a re-implementation of NASA GRC's EZ FSP Sizer [mason2025trades].
EZ FSP Sizer combines first-order physics with hardware-based curve fits and
mass models. TLECUIL contains only first-principles physics plus inputs you
supply. It has no mass, shield, dose, or transient model (Theory Sec. 14 [sec:limits]).

Scope of the physics
--------------------
  1. Electrical chain: bus <- PMAD <- parasitics <- converter   (Theory Sec. 3 [sec:power])
  2. Converter: first law + fixed efficiency or fraction of Carnot  (Sec. 4 [sec:eta])
  3. Reactor-to-converter heat losses                            (Sec. 4.4 [sec:loss])
  4. Main radiator: isothermal or pumped fluid loop, with panel
     effectiveness, two-sided panels and area margin             (Sec. 5 [sec:rad])
  5. Electronics radiator                                         (Sec. 6 [sec:elec])
  6. U-235 fissioned (energy-equivalent)                          (Sec. 7 [sec:fuel])
  7. Cold-end temperature that minimises radiator area            (Sec. 8 [sec:opt])

Units: kW (kWe for electric, kWt for thermal), K, m^2, kg, years.

Usage
-----
    python tlecuil.py                       # default 40 kWe case
    python tlecuil.py --power 100           # any power; no upper limit
    python tlecuil.py --power 100 --plots   # trade-study PNGs
    python tlecuil.py --export-json cases.json  # cases for the HTML parity check
    python benchmarks.py                    # compare against NASA design points
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict, replace
from typing import Literal

__version__ = "0.2.0"

# ---------------------------------------------------------------------------
# 0. PHYSICAL CONSTANTS -- CODATA 2018 [codata2018]
# ---------------------------------------------------------------------------
# sigma, N_A and the MeV->J factor are exact or recommended values from the
# 2018 CODATA adjustment. The U-235 atomic mass is the standard tabulated
# value. Keeping them here means one place to audit.

STEFAN_BOLTZMANN = 5.670374419e-8    # sigma, W m^-2 K^-4      [codata2018]
AVOGADRO = 6.02214076e23             # N_A, mol^-1 (exact)      [codata2018]
MEV_TO_J = 1.602176634e-13           # 1 MeV in J (exact e)     [codata2018]
M_U235 = 235.0439299                 # g/mol, U-235 atomic mass
SECONDS_PER_YEAR = 365.25 * 86400.0  # Julian year, s


# ---------------------------------------------------------------------------
# 1. INPUTS
# ---------------------------------------------------------------------------
# A frozen dataclass: every input is named, typed and has a default, and an
# input set cannot be changed by accident. `replace(inp, field=value)` makes
# a modified copy, which is how trade studies and benchmarks work.
#
# DEFAULTS: chosen to sit inside the range of published NASA FSP studies,
# with the source noted. They describe a generic 40 kWe lunar system, not
# any specific NASA design. Benchmarks (benchmarks.py) set inputs that
# reproduce specific published designs.

@dataclass(frozen=True)
class FSPInputs:
    # --- Electrical chain (Theory Sec. 3 [sec:power]) ----------------------------------
    delivered_power_kwe: float = 40.0
    # 40 kWe is the NASA FSP net requirement in [fsp2010icd] and
    # [mason2025trades]. NASA's 2025 solicitation moved to >= 100 kWe
    # [nasa2025fsp]; there is no upper limit in this code.
    pmad_efficiency: float = 0.90
    # Bus power / PMAD input. ~90% (Brayton) and 83% (Stirling) in
    # [mason2025trades]; 87% in [oleson2022deployable].
    parasitic_fraction: float = 0.0
    # Share of converter output used by loads that scale with power
    # (e.g. pumps). The 2010 reference had ~5 kWe of parasitics on
    # 48 kWe gross, about 10% [fsp2010icd].
    fixed_parasitic_kwe: float = 0.0
    # Loads that do NOT scale with power (housekeeping, instruments).
    # [mason2025trades] allows up to 2 kWe of housekeeping power.

    # --- Power conversion (Theory Sec. 4 [sec:eta]) ----------------------------------
    efficiency_mode: Literal["fixed", "carnot"] = "carnot"
    conversion_efficiency: float = 0.26
    # Used when mode == "fixed". 26% projected for the 12 kWe Stirling
    # [fsp2010icd]; 21.7-23.3% measured in TDU tests [briggs2016tdu].
    fraction_of_carnot: float = 0.55
    # Used when mode == "carnot". Values implied by published designs
    # range from about 0.3 (Brayton) to 0.57 (Stirling), depending on
    # how Th and Tc are defined (Theory Table 11.4 [tab:carnot]).
    hot_temperature_k: float = 850.0
    cold_temperature_k: float = 425.0
    # Effective cycle temperatures. 778 K / 425 K in [fsp2010icd];
    # 850 K is the NaK reactor outlet in the same design.
    reactor_heat_loss_fraction: float = 0.0
    # Share of reactor heat that never reaches the converter (leakage
    # from reactor and piping). 18% reactor losses in the heat-pipe
    # design of [oleson2022deployable]. This heat is not sized here.

    # --- Main radiator (Theory Sec. 5 [sec:rad]) -------------------------------------
    hx_delta_t_k: float = 5.0
    # Drop from converter cold end to radiator coolant inlet. 425 K cold
    # end vs 420 K water outlet in [fsp2010icd].
    radiator_delta_t_k: float = 30.0
    # Inlet-to-outlet coolant drop: 420 -> 390 K water in [fsp2010icd].
    # Set to 0 for an isothermal (ideal heat-pipe) radiator.
    sink_temperature_k: float = 250.0
    # Effective sink. 250 K used for the 2010 design [fsp2010icd]; 270 K
    # SOW value [mason2025trades]. Real sinks vary from ~100 K to ~320 K
    # with site and time of day [mason2025trades].
    emissivity: float = 0.85
    # End-of-life emissivity used in [mason2025trades].
    radiating_sides: int = 2
    # NASA FSP panels are two-sided vertical panels [fsp2010icd].
    panel_effectiveness: float = 0.76
    # Lumps fin efficiency, heat-pipe and facesheet temperature drops,
    # and non-ideal view factors. 1.0 is an ideal panel. 0.76 is
    # back-calculated from the 2010 NASA design alone (0.756) and then
    # tested against independent designs (Theory Sec. 11 [sec:bench], benchmarks.py).
    # It is a calibration factor, not a measured property.
    area_margin: float = 0.10
    # NASA practice: 10% area margin [fsp2010icd; mason2025trades].

    # --- Electronics radiator (Theory Sec. 6 [sec:elec]) ------------------------------
    electronics_radiator_temperature_k: float = 333.0
    # Electronics radiator temperature in [mason2025trades].

    # --- Fuel (Theory Sec. 7 [sec:fuel]) ----------------------------------------------
    full_power_years: float = 10.0
    # 10-year requirement in [mason2025trades]; 8 years in [fsp2010icd].
    energy_per_fission_mev: float = 200.0
    # Recoverable energy per U-235 fission, ~200 MeV [lamarsh2017].


# ---------------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------------
# A formula does not crash on non-physical inputs; it returns a wrong
# number. So each physical constraint is checked first and explained.

def validate(inp: FSPInputs) -> None:
    def need(ok: bool, msg: str) -> None:
        if not ok:
            raise ValueError(msg)

    need(inp.delivered_power_kwe > 0, "Delivered power must be positive.")
    need(0 < inp.pmad_efficiency <= 1, "PMAD efficiency must be in (0, 1].")
    need(0 <= inp.parasitic_fraction < 1, "Parasitic fraction must be in [0, 1).")
    need(inp.fixed_parasitic_kwe >= 0, "Fixed parasitic load cannot be negative.")
    need(0 <= inp.reactor_heat_loss_fraction < 1, "Reactor heat-loss fraction must be in [0, 1).")
    need(inp.efficiency_mode in ("fixed", "carnot"), "Mode must be 'fixed' or 'carnot'.")
    need(inp.hot_temperature_k > inp.cold_temperature_k > 0,
         "Need hot temperature > cold temperature > 0 K.")
    carnot = 1.0 - inp.cold_temperature_k / inp.hot_temperature_k
    if inp.efficiency_mode == "fixed":
        need(0 < inp.conversion_efficiency < 1, "Conversion efficiency must be in (0, 1).")
        # Second law (Theory eq. 4.3 [eq:carnot]): no engine beats Carnot.
        need(inp.conversion_efficiency <= carnot,
             f"Fixed efficiency {inp.conversion_efficiency:.3f} exceeds the Carnot "
             f"limit {carnot:.3f} for these temperatures (second law).")
    else:
        need(0 < inp.fraction_of_carnot <= 1, "Fraction of Carnot must be in (0, 1].")
    need(0 < inp.emissivity <= 1, "Emissivity must be in (0, 1].")
    need(inp.radiating_sides in (1, 2), "Radiating sides must be 1 or 2.")
    need(0 < inp.panel_effectiveness <= 1, "Panel effectiveness must be in (0, 1].")
    need(inp.area_margin >= 0, "Area margin cannot be negative.")
    need(inp.hx_delta_t_k >= 0 and inp.radiator_delta_t_k >= 0,
         "Temperature drops cannot be negative.")
    need(inp.sink_temperature_k >= 0, "Sink temperature cannot be negative.")
    t_out = inp.cold_temperature_k - inp.hx_delta_t_k - inp.radiator_delta_t_k
    need(t_out > inp.sink_temperature_k,
         "Radiator outlet must stay hotter than the sink, or it cannot reject heat.")
    need(inp.electronics_radiator_temperature_k > inp.sink_temperature_k,
         "Electronics radiator must be hotter than the sink.")
    need(inp.full_power_years >= 0, "Mission life cannot be negative.")
    need(inp.energy_per_fission_mev > 0, "Energy per fission must be positive.")


# ---------------------------------------------------------------------------
# 3. CONVERSION EFFICIENCY (Theory Sec. 4 [sec:eta])
# ---------------------------------------------------------------------------
# Carnot limit (Theory eq. 4.3 [eq:carnot]): eta_C = 1 - Tc/Th, the second-law ceiling
# for any heat engine [smith2001thermo]. Real Stirling and Brayton machines
# reach a fraction f of it (Theory eq. 4.4 [eq:fc]):  eta = f (1 - Tc/Th).
# This links the radiator to the reactor: a hotter cold end shrinks the
# radiator (T^4) but lowers eta, so the reactor must make more heat.

def conversion_efficiency(inp: FSPInputs) -> float:
    if inp.efficiency_mode == "fixed":
        return inp.conversion_efficiency
    return inp.fraction_of_carnot * (1.0 - inp.cold_temperature_k / inp.hot_temperature_k)


# ---------------------------------------------------------------------------
# 4. RADIATOR AREA (Theory Sec. 5 [sec:rad])
# ---------------------------------------------------------------------------
# Net flux leaving a gray, diffuse surface at T toward an effective sink at
# Ts (Theory eq. 5.1 [eq:flux]) [howell2021radiation; gilmore2002thermal]:
#     q'' = eps * sigma * (T^4 - Ts^4)              [W per m^2 of surface]
#
# This is per square metre of RADIATING SURFACE. A two-sided panel has two
# radiating faces, so its planform (the panel as built) is half the
# radiating area. NASA quotes the total of both faces [fsp2010icd;
# mason2025trades], so TLECUIL reports that as the primary result.
#
# (a) Isothermal surface (Theory eq. 5.5 [eq:aiso]):
#     A = Q / (eta_p * eps * sigma * (T^4 - Ts^4))
# (b) Pumped fluid loop, coolant cooling from T_in to T_out (Theory
#     eqs. 5.6-5.7 [eq:strip, eq:afluid]). An energy balance on a strip dA gives
#     -mdot cp dT = eta_p eps sigma (T^4 - Ts^4) dA, which integrates to
#     A = mdot cp / (eta_p eps sigma) * Integral_{T_out}^{T_in} dT / (T^4 - Ts^4),
#     with mdot cp = Q / (T_in - T_out).
# The integral has a closed form (Theory eq. 5.8 [eq:antider]) and, for small sinks, a
# series form (eq. 5.9 [eq:series]) that avoids round-off; both are coded below.
#
# eta_p is the panel effectiveness (Theory Sec. 5.3 [sec:etap]). A fin-efficiency
# factor multiplying the radiated flux is standard in space-radiator
# models [wang2021entropy; gilmore2002thermal].

def _antiderivative(t: float, a: float) -> float:
    """F(T) with dF/dT = 1/(T^4 - a^4), valid for T > a >= 0 (Theory eq. 5.8 [eq:antider])."""
    if a < 1e-9:                                   # sink at ~0 K: integral of T^-4
        return -1.0 / (3.0 * t ** 3)
    return (math.log((t - a) / (t + a)) / (4.0 * a ** 3)
            - math.atan(t / a) / (2.0 * a ** 3))


def _integral_inv_t4(t_in: float, t_out: float, a: float) -> float:
    """Integral of dT / (T^4 - a^4) from t_out to t_in (Theory eqs. 5.8-5.9 [eq:antider, eq:series]).

    The closed form (eq. 5.8 [eq:antider]) subtracts two terms of size ~1/a^3. When the
    sink a is small compared with T they nearly cancel and round-off destroys
    the answer (for a = 1e-3 K it returns a negative area). So for a/T_out <
    0.5 we use the convergent series 1/(T^4 - a^4) = sum_k a^(4k) / T^(4k+4),
    integrated term by term (eq. 5.9 [eq:series]). Each term is at most 0.5^4 = 1/16 of
    the previous one, so 20 terms reach machine precision. The two methods
    agree to ~1e-12 where both are valid; tests check both against SciPy.
    """
    if a < 0.5 * t_out:
        total, x = 0.0, 1.0                        # x = a^(4k)
        for k in range(20):
            n = 4 * k + 3
            total += x / n * (t_out ** -n - t_in ** -n)
            x *= a ** 4
        return total
    return _antiderivative(t_in, a) - _antiderivative(t_out, a)


def radiating_area_m2(q_watts: float, t_in: float, t_out: float, t_sink: float,
                      emissivity: float, panel_effectiveness: float) -> float:
    """Total radiating surface (all faces) needed to reject q_watts, before margin."""
    if q_watts <= 0:
        return 0.0
    k = panel_effectiveness * emissivity * STEFAN_BOLTZMANN   # W m^-2 K^-4
    if t_in - t_out < 1e-6:                                   # case (a)
        return q_watts / (k * (t_in ** 4 - t_sink ** 4))
    mdot_cp = q_watts / (t_in - t_out)                        # case (b), W/K
    return mdot_cp * _integral_inv_t4(t_in, t_out, t_sink) / k


# ---------------------------------------------------------------------------
# 5. THE SIZING CALCULATION
# ---------------------------------------------------------------------------
# Walk the energy chain backwards from the user requirement:
#   bus <- PMAD <- (minus parasitics) <- converter <- (minus losses) <- reactor
# and record where every lost kilowatt goes as heat.

@dataclass(frozen=True)
class FSPResults:
    conversion_efficiency: float
    carnot_efficiency: float
    converter_output_kwe: float          # gross electric at converter terminals
    parasitic_load_kwe: float            # proportional + fixed
    pmad_loss_kw: float
    converter_heat_input_kwt: float      # heat delivered to the converter
    reactor_thermal_kwt: float           # heat produced by fission
    reactor_heat_loss_kwt: float         # not rejected by the main radiator
    converter_heat_rejection_kwt: float  # main radiator load
    electronics_heat_kwt: float          # PMAD loss + parasitics
    system_efficiency: float             # delivered / reactor thermal
    radiator_inlet_k: float
    radiator_outlet_k: float
    radiator_area_m2: float              # total radiating area, incl. margin
    radiator_planform_m2: float          # panel size = area / sides
    electronics_radiator_area_m2: float  # incl. margin
    u235_fissioned_kg: float
    energy_balance_error_kw: float       # self-check, should be ~0


def size_fsp(inp: FSPInputs) -> FSPResults:
    validate(inp)

    # Theory eq. 3.2 [eq:pconv]: P_del = eta_PMAD * [(1 - f_par) P_conv - P_fix]
    #   => P_conv = (P_del / eta_PMAD + P_fix) / (1 - f_par)
    p_conv = ((inp.delivered_power_kwe / inp.pmad_efficiency + inp.fixed_parasitic_kwe)
              / (1.0 - inp.parasitic_fraction))
    p_par = inp.parasitic_fraction * p_conv + inp.fixed_parasitic_kwe
    p_pmad_loss = (p_conv - p_par) - inp.delivered_power_kwe            # Theory eq. 3.3 [eq:lpmad]

    # Theory eqs. 4.1-4.2 [eq:firstlaw, eq:qin]: first law around the converter
    eta = conversion_efficiency(inp)
    q_in = p_conv / eta
    q_conv = q_in - p_conv

    # Theory eq. 4.5 [eq:loss]: reactor heat that leaks before reaching the converter
    p_th = q_in / (1.0 - inp.reactor_heat_loss_fraction)
    q_loss = p_th - q_in

    # Theory Sec. 6 [sec:elec]: electronics heat. Assumption: all parasitics end as
    # low-temperature heat. This is conservative for TOTAL radiator area (the
    # electronics radiator is colder) but slightly under-loads the main
    # radiator if, in reality, pump work heats the main coolant.
    q_elec = p_pmad_loss + p_par

    # Theory Sec. 5 [sec:rad]: main radiator
    t_in = inp.cold_temperature_k - inp.hx_delta_t_k                   # Theory eq. 5.4 [eq:tin]
    t_out = t_in - inp.radiator_delta_t_k
    margin = 1.0 + inp.area_margin
    a_main = margin * radiating_area_m2(q_conv * 1e3, t_in, t_out, inp.sink_temperature_k,
                                        inp.emissivity, inp.panel_effectiveness)
    t_e = inp.electronics_radiator_temperature_k
    a_elec = margin * radiating_area_m2(q_elec * 1e3, t_e, t_e, inp.sink_temperature_k,
                                        inp.emissivity, inp.panel_effectiveness)

    # Theory eq. 7.1 [eq:mf]: fissions needed for P_th over the mission
    fission_rate = p_th * 1e3 / (inp.energy_per_fission_mev * MEV_TO_J)   # s^-1
    n_fissions = fission_rate * inp.full_power_years * SECONDS_PER_YEAR
    m_fissioned_kg = n_fissions / AVOGADRO * M_U235 / 1e3

    # Theory eq. 6.2 [eq:balance]: everything the reactor makes leaves as bus power or heat
    balance = p_th - (inp.delivered_power_kwe + q_conv + q_elec + q_loss)

    return FSPResults(
        conversion_efficiency=eta,
        carnot_efficiency=1.0 - inp.cold_temperature_k / inp.hot_temperature_k,
        converter_output_kwe=p_conv,
        parasitic_load_kwe=p_par,
        pmad_loss_kw=p_pmad_loss,
        converter_heat_input_kwt=q_in,
        reactor_thermal_kwt=p_th,
        reactor_heat_loss_kwt=q_loss,
        converter_heat_rejection_kwt=q_conv,
        electronics_heat_kwt=q_elec,
        system_efficiency=inp.delivered_power_kwe / p_th,
        radiator_inlet_k=t_in,
        radiator_outlet_k=t_out,
        radiator_area_m2=a_main,
        radiator_planform_m2=a_main / inp.radiating_sides,
        electronics_radiator_area_m2=a_elec,
        u235_fissioned_kg=m_fissioned_kg,
        energy_balance_error_kw=balance,
    )


# ---------------------------------------------------------------------------
# 6. COLD-END TEMPERATURE THAT MINIMISES RADIATOR AREA (Theory Sec. 8 [sec:opt])
# ---------------------------------------------------------------------------
# Only meaningful in "carnot" mode, where eta depends on Tc. Golden-section
# search between "just above the sink" and "just below Th". The same
# algorithm is written in the HTML page so both give identical answers.
# It is tested against the closed form (Theory eq. 8.2 [eq:tauopt]), which reduces to
# the classic Tc/Th = 3/4 for an ideal engine [smith2001thermo].
#
# Minimum AREA is not minimum MASS: at the area optimum the efficiency is
# low and the reactor is large. Treat it as one side of the trade.

def optimal_cold_temperature(inp: FSPInputs, tol_k: float = 1e-3) -> tuple[float, float]:
    if inp.efficiency_mode != "carnot":
        raise ValueError("Optimum Tc is only defined when efficiency depends on Tc.")
    lo = inp.sink_temperature_k + inp.hx_delta_t_k + inp.radiator_delta_t_k + 1.0
    hi = inp.hot_temperature_k - 1.0

    def area(tc: float) -> float:
        return size_fsp(replace(inp, cold_temperature_k=tc)).radiator_area_m2

    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = area(c), area(d)
    while b - a > tol_k:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = area(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = area(d)
    tc = 0.5 * (a + b)
    return tc, area(tc)


def analytic_optimum_ratio(f: float) -> float:
    """Tc/Th that minimises area for Ts = 0, no drops, isothermal (Theory eq. 8.2 [eq:tauopt])."""
    b = 5.0 - 8.0 * f
    return (-b + math.sqrt(b * b + 64.0 * f * (1.0 - f))) / (8.0 * f)


# ---------------------------------------------------------------------------
# 7. TRADE STUDIES
# ---------------------------------------------------------------------------

def sweep(inp: FSPInputs, field: str, values) -> list[tuple[float, FSPResults]]:
    """Re-run size_fsp for each value of one input; skip non-physical points."""
    out = []
    for v in values:
        try:
            out.append((v, size_fsp(replace(inp, **{field: v}))))
        except ValueError:
            pass
    return out


def make_plots(inp: FSPInputs, prefix: str = "tlecuil") -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    files = []
    lo = inp.sink_temperature_k + inp.hx_delta_t_k + inp.radiator_delta_t_k + 5
    pts = sweep(inp, "cold_temperature_k", np.linspace(lo, inp.hot_temperature_k - 5, 300))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot([p[0] for p in pts], [p[1].radiator_area_m2 for p in pts], lw=2)
    axes[0].set(xlabel="Converter cold-end temperature $T_c$ (K)",
                ylabel="Radiating area, incl. margin (m$^2$)",
                title=f"Radiator area, {inp.delivered_power_kwe:g} kWe", yscale="log")
    axes[1].plot([p[0] for p in pts], [p[1].reactor_thermal_kwt for p in pts], lw=2)
    axes[1].set(xlabel="Converter cold-end temperature $T_c$ (K)",
                ylabel="Reactor thermal power (kWt)", title="Reactor power", yscale="log")
    if inp.efficiency_mode == "carnot":
        tc, a = optimal_cold_temperature(inp)
        axes[0].axvline(tc, ls="--", c="gray")
        axes[0].annotate(f"min {a:.0f} m$^2$ at {tc:.0f} K", (tc, a),
                         xytext=(10, 30), textcoords="offset points")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    files.append(f"{prefix}_tc_trade.png")
    fig.savefig(files[-1], dpi=150)
    plt.close(fig)
    return files


# ---------------------------------------------------------------------------
# 8. REPORTING AND COMMAND LINE
# ---------------------------------------------------------------------------

def report(inp: FSPInputs, r: FSPResults) -> str:
    rows = [
        ("Delivered power", f"{inp.delivered_power_kwe:.2f} kWe"),
        ("Conversion efficiency", f"{r.conversion_efficiency:.1%} (Carnot {r.carnot_efficiency:.1%})"),
        ("Converter output (gross)", f"{r.converter_output_kwe:.2f} kWe"),
        ("Parasitic load", f"{r.parasitic_load_kwe:.2f} kWe"),
        ("Reactor thermal power", f"{r.reactor_thermal_kwt:.1f} kWt"),
        ("Main radiator heat load", f"{r.converter_heat_rejection_kwt:.1f} kWt"),
        ("Electronics heat", f"{r.electronics_heat_kwt:.2f} kWt"),
        ("Reactor/piping heat loss", f"{r.reactor_heat_loss_kwt:.2f} kWt"),
        ("System efficiency", f"{r.system_efficiency:.1%}"),
        ("Radiator inlet / outlet", f"{r.radiator_inlet_k:.0f} / {r.radiator_outlet_k:.0f} K"),
        ("Main radiator, radiating area", f"{r.radiator_area_m2:.1f} m^2 (incl. margin)"),
        ("Main radiator, panel planform", f"{r.radiator_planform_m2:.1f} m^2"),
        ("Electronics radiator area", f"{r.electronics_radiator_area_m2:.1f} m^2"),
        ("U-235 fissioned over life", f"{r.u235_fissioned_kg:.3f} kg"),
        ("Energy balance error", f"{r.energy_balance_error_kw:.2e} kW"),
    ]
    w = max(len(k) for k, _ in rows)
    return "\n".join(f"{k:<{w}}  {v}" for k, v in rows)


PARITY_CASES = {
    "default_40kWe": FSPInputs(),
    "100kWe": FSPInputs(delivered_power_kwe=100.0),
    "fixed_eta_isothermal": FSPInputs(
        delivered_power_kwe=75, efficiency_mode="fixed", conversion_efficiency=0.27,
        cold_temperature_k=450, hx_delta_t_k=0, radiator_delta_t_k=0,
        sink_temperature_k=210, radiating_sides=1, panel_effectiveness=1.0, area_margin=0.0),
    "losses_and_fixed_loads": FSPInputs(
        delivered_power_kwe=60, parasitic_fraction=0.05, fixed_parasitic_kwe=2.0,
        reactor_heat_loss_fraction=0.1, radiator_delta_t_k=80, hot_temperature_k=1000),
    "1MWe": FSPInputs(delivered_power_kwe=1000.0, hot_temperature_k=1150.0),
    "zero_sink": FSPInputs(sink_temperature_k=0.0),
    "tiny_sink": FSPInputs(sink_temperature_k=1e-3),
    "fixed_eta_fluid_loop": FSPInputs(efficiency_mode="fixed", conversion_efficiency=0.2,
                                      radiator_delta_t_k=100, sink_temperature_k=200),
    "brayton_100kWe_preset": FSPInputs(
        delivered_power_kwe=100, pmad_efficiency=0.90, fixed_parasitic_kwe=2,
        efficiency_mode="fixed", conversion_efficiency=0.219, hot_temperature_k=1100,
        cold_temperature_k=504, hx_delta_t_k=0, radiator_delta_t_k=125, sink_temperature_k=270),
}


def main() -> None:
    ap = argparse.ArgumentParser(description="TLECUIL first-order FSP sizer")
    ap.add_argument("--power", type=float, help="delivered power, kWe (no upper limit)")
    ap.add_argument("--mode", choices=["fixed", "carnot"])
    ap.add_argument("--eta", type=float, help="fixed conversion efficiency")
    ap.add_argument("--tc", type=float, help="cold-end temperature, K")
    ap.add_argument("--th", type=float, help="hot-end temperature, K")
    ap.add_argument("--sink", type=float, help="effective sink temperature, K")
    ap.add_argument("--plots", action="store_true", help="write trade-study PNGs")
    ap.add_argument("--export-json", metavar="FILE", help="write cases for the HTML parity check")
    args = ap.parse_args()

    if args.export_json:
        data = {}
        for name, case in PARITY_CASES.items():
            entry = {"inputs": asdict(case), "results": asdict(size_fsp(case))}
            if case.efficiency_mode == "carnot":
                entry["optimum_tc_k"], entry["optimum_area_m2"] = optimal_cold_temperature(case)
            data[name] = entry
        with open(args.export_json, "w") as fh:
            json.dump(data, fh, indent=2)
        print(f"wrote {args.export_json}")
        return

    overrides = {k: v for k, v in {
        "delivered_power_kwe": args.power, "efficiency_mode": args.mode,
        "conversion_efficiency": args.eta, "cold_temperature_k": args.tc,
        "hot_temperature_k": args.th, "sink_temperature_k": args.sink}.items() if v is not None}
    if args.eta is not None and args.mode is None:
        overrides["efficiency_mode"] = "fixed"
    inp = replace(FSPInputs(), **overrides)
    res = size_fsp(inp)
    print(f"TLECUIL {__version__}\n")
    print(report(inp, res))
    if inp.efficiency_mode == "carnot":
        tc, a = optimal_cold_temperature(inp)
        print(f"\nMinimum-area cold-end temperature: {tc:.1f} K ({a:.1f} m^2)")
    if args.plots:
        for f in make_plots(inp):
            print("wrote", f)


if __name__ == "__main__":
    main()
