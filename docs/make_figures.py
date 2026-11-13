"""
make_figures.py -- generate the figures used in tlecuil_theory.tex.

Run from the docs/ folder:   python make_figures.py
Every figure is computed by tlecuil.py itself, so the pictures in the paper
always match what the code does. Re-run after any change to the code.
"""
import os
import sys
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from tlecuil import (FSPInputs, size_fsp, sweep, optimal_cold_temperature,   # noqa: E402
                     STEFAN_BOLTZMANN, _integral_inv_t4)
import benchmarks                                                            # noqa: E402

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
BLUE, ORANGE, GREEN, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#7c7a74"


def fig_flux_vs_temperature():
    """Why temperature matters so much: radiated flux grows as T^4."""
    t = np.linspace(250, 1000, 300)
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    for ts, c, lab in ((0, BLUE, "sink 0 K (deep space)"), (250, ORANGE, "sink 250 K (lunar)")):
        q = 0.85 * STEFAN_BOLTZMANN * (t ** 4 - ts ** 4) / 1e3
        ax.plot(t, np.where(t > ts, q, np.nan), color=c, lw=2, label=lab)
    for tt in (400, 800):
        q = 0.85 * STEFAN_BOLTZMANN * tt ** 4 / 1e3
        ax.plot(tt, q, "o", color=BLUE)
        ax.annotate(f"{tt} K: {q:.1f} kW/m$^2$", (tt, q), xytext=(10, -12), textcoords="offset points")
    ax.set(xlabel="Radiator surface temperature $T$ (K)",
           ylabel="Net flux per m$^2$ (kW/m$^2$)", title=r"$q''=\varepsilon\sigma(T^4-T_s^4)$, $\varepsilon=0.85$")
    ax.legend(frameon=False)
    fig.savefig(os.path.join(OUT, "flux_vs_T.pdf"))
    plt.close(fig)


def fig_tc_trade():
    """The central trade: hotter cold end -> smaller radiator, bigger reactor."""
    inp = FSPInputs()
    lo = inp.sink_temperature_k + inp.hx_delta_t_k + inp.radiator_delta_t_k + 5
    pts = sweep(inp, "cold_temperature_k", np.linspace(lo, inp.hot_temperature_k - 5, 400))
    tc = [p[0] for p in pts]
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
    ax[0].plot(tc, [p[1].radiator_area_m2 for p in pts], color=BLUE, lw=2)
    ax[1].plot(tc, [p[1].reactor_thermal_kwt for p in pts], color=ORANGE, lw=2)
    t_opt, a_opt = optimal_cold_temperature(inp)
    ax[0].axvline(t_opt, ls="--", color=GRAY)
    ax[0].annotate(f"minimum area\n{a_opt:.0f} m$^2$ at {t_opt:.0f} K", (t_opt, a_opt),
                   xytext=(-110, 60), textcoords="offset points", arrowprops=dict(arrowstyle="->"))
    r = size_fsp(inp)
    for a, y in ((ax[0], r.radiator_area_m2), (ax[1], r.reactor_thermal_kwt)):
        a.plot(inp.cold_temperature_k, y, "o", color="k")
        a.annotate("default design", (inp.cold_temperature_k, y), xytext=(8, 8), textcoords="offset points")
        a.set_yscale("log")
        a.set_xlabel("Cold-end temperature $T_c$ (K)")
    ax[0].set(ylabel="Radiating area (m$^2$)", title="Radiator")
    ax[1].set(ylabel="Reactor thermal power (kWt)", title="Reactor")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "tc_trade.pdf"))
    plt.close(fig)


def fig_fluid_loop_profile():
    """Coolant temperature along a pumped-loop radiator (from Theory eq. 5.7 [eq:afluid])."""
    t_in, t_out, ts = 504.0, 379.0, 270.0      # the 2025 Brayton loop (Theory Sec. 11 [sec:bench])
    t = np.linspace(t_in, t_out, 200)
    frac = np.array([_integral_inv_t4(t_in, x, ts) if x < t_in else 0.0 for x in t])
    frac = frac / frac[-1]
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(frac, t, color=BLUE, lw=2, label="pumped loop (exact)")
    ax.plot([0, 1], [t_in, t_out], ls="--", color=GRAY, label="straight line (for comparison)")
    ax.set(xlabel="Fraction of radiator area from the inlet", ylabel="Coolant temperature (K)",
           title="Coolant cools fastest near the hot inlet")
    ax.legend(frameon=False)
    fig.savefig(os.path.join(OUT, "fluid_loop_profile.pdf"))
    plt.close(fig)


def fig_benchmark_parity():
    """Predicted vs published radiator areas (ideal panel and calibrated)."""
    res = benchmarks.run()
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    pub = [r["area_pub"] for r in res["radiator"]]
    ax.plot([100, 200], [100, 200], color=GRAY, lw=1)
    ax.fill_between([100, 200], [85, 170], [115, 230], color=GRAY, alpha=0.12, label="$\\pm$15% band")
    ax.plot(pub, [r["area_ideal"] for r in res["radiator"]], "s", color=ORANGE, label="ideal panel ($\\eta_p=1$)")
    ax.plot(pub, [r["area_pred"] for r in res["radiator"]], "o", color=BLUE,
            label=f"calibrated ($\\eta_p={res['eta_p_calibrated']:.3f}$)")
    offs = {"B1": (6, -10), "B3": (-18, 6)}
    for r in res["radiator"]:
        ax.annotate(r["key"], (r["area_pub"], r["area_pred"]), xytext=offs.get(r["key"], (5, -3)),
                    textcoords="offset points")
    ax.set(xlim=(100, 200), ylim=(100, 200), xlabel="Published area (m$^2$)",
           ylabel="TLECUIL area (m$^2$)", title="Radiator benchmark")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(os.path.join(OUT, "benchmark_parity.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_flux_vs_temperature()
    fig_tc_trade()
    fig_fluid_loop_profile()
    fig_benchmark_parity()
    print("figures written to", OUT)
