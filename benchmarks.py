"""
benchmarks.py -- compare TLECUIL against published NASA FSP design points
==========================================================================

Run:  python benchmarks.py            (prints tables, writes JSON + LaTeX)

What this file does, and why it is set up this way
--------------------------------------------------
A benchmark is only meaningful if (1) each input is taken from the source,
or the assumption is written down; (2) the comparison is like for like (the
same area definition, margin and heat load); (3) any value tuned to match a
source is NOT then counted as agreement with that same source.

So the radiator comparison is done in three steps:

  Step 1 -- IDEAL: an ideal panel (effectiveness = 1), every other input
            from the source. This shows the raw first-principles gap.
  Step 2 -- CALIBRATE on ONE design only (B1, the 2010 NASA reference) to
            find the single panel effectiveness eta_p that reproduces its
            published area.
  Step 3 -- PREDICT the other, independent designs (B2-B4) with that eta_p
            fixed. Agreement there is real evidence. Agreement on B1 is not:
            it was fitted.

The radiator check uses the PUBLISHED heat load, so it tests only the
radiator physics. The system checks (gross power, reactor power) test the
electrical chain and the first law separately.

Sources (keys match docs/references.bib):
  fsp2010icd        NASA/TM-2010-216772, FSP Initial Concept Definition
  mason2011arch     Mason & Poston, NASA/TM-2011-216819
  mason2025trades   Mason, Kaldon, Corbisiero, Rao, NETS 2025
  oleson2022deployable  Oleson et al., Deployable 40 kWe Lunar FSP Concept
  briggs2016tdu     Briggs et al., FSP TDU test results (NASA/TM-2016-219382)
  dhar1999sse       Dhar, Stirling Space Engine Program final report
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, replace

from tlecuil import FSPInputs, size_fsp, radiating_area_m2

TOL_RADIATOR = 0.15     # +/-15 % acceptance for first-order radiator predictions
TOL_POWER = 0.05        # +/-5 % for electrical-chain and first-law quantities


# ---------------------------------------------------------------------------
# Radiator benchmark cases. Every number is from the cited source unless
# listed under "assumed". Areas are the TOTAL radiating area (both faces)
# including the source's margin.
# ---------------------------------------------------------------------------

@dataclass
class RadiatorCase:
    key: str
    label: str
    source: str
    q_kwt: float                  # published main-radiator heat load
    t_in: float                   # coolant inlet, K
    t_out: float                  # coolant outlet, K
    t_sink: float                 # effective sink, K
    emissivity: float
    margin: float
    area_pub_m2: float            # published total radiating area
    assumed: list = field(default_factory=list)


RADIATOR_CASES = [
    RadiatorCase(
        key="B1", label="2010 FSP reference, Stirling, 40 kWe",
        source="fsp2010icd; mason2011arch",
        q_kwt=140.0, t_in=420.0, t_out=390.0, t_sink=250.0,
        emissivity=0.85, margin=0.10, area_pub_m2=185.0,
        assumed=["emissivity 0.85 (not given in the source; product eta_p*eps is what is calibrated)",
                 "heat load is '~140 kWt' in the source"]),
    RadiatorCase(
        key="B2", label="2025 GRD, heat-pipe reactor + Stirling, 40 kWe",
        source="mason2025trades",
        q_kwt=119.0, t_in=433.0, t_out=403.0, t_sink=270.0,
        emissivity=0.85, margin=0.10, area_pub_m2=140.0,
        assumed=["'radiator operating range 403-433 K' read as loop inlet/outlet",
                 "sink 270 K = SOW value; the source does not say which sink sized the radiator",
                 "area = 14 panels x 10 m^2 (both faces)"]),
    RadiatorCase(
        key="B3", label="2025 GRD, gas-cooled reactor + Brayton, 40 kWe",
        source="mason2025trades",
        q_kwt=180.0, t_in=504.0, t_out=379.0, t_sink=270.0,
        emissivity=0.85, margin=0.10, area_pub_m2=180.0,
        assumed=["'operating between 379 and 504 K' read as loop inlet/outlet (pumped water loop)",
                 "sink 270 K = SOW value",
                 "area = 18 panels x 10 m^2 (both faces)"]),
    RadiatorCase(
        key="B4", label="Deployable 40 kWe concept, Stirling, south pole",
        source="oleson2022deployable",
        q_kwt=126.4, t_in=420.0, t_out=370.0,
        t_sink=round((0.5 * 220.0 ** 4) ** 0.25, 1), emissivity=0.84, margin=0.0,
        area_pub_m2=133.4,
        assumed=["sink built from the source's view factors: 0.5 to a 220 K surface, 0.5 to deep space",
                 "margin not stated; 0 assumed",
                 "source does not say whether 133.4 m^2 is planform or total radiating area"]),
]


def radiator_area(case: RadiatorCase, eta_p: float) -> float:
    return (1 + case.margin) * radiating_area_m2(
        case.q_kwt * 1e3, case.t_in, case.t_out, case.t_sink, case.emissivity, eta_p)


def calibrate_eta_p(case: RadiatorCase) -> float:
    """eta_p that reproduces the published area exactly (area is proportional to 1/eta_p)."""
    return radiator_area(case, 1.0) / case.area_pub_m2


# ---------------------------------------------------------------------------
# System-level cases: electrical chain and first law
# ---------------------------------------------------------------------------

def system_cases(p_fix_2025=2.0):
    """(key, label, source, inputs, {metric: published}, identities)

    S1 gross output is an identity (40 + 3 + 5 = 48 by construction), so it is
    reported but not counted as agreement. For S2/S3 the 2 kWe fixed load is an
    INTERPRETATION: [mason2025trades] allows "up to 2 kWe" of housekeeping,
    supplied in part by an auxiliary source; the source does not say it is
    drawn from the converter. Results are reported with and without it.
    """
    return [
        ("S1", "2010 reference: 3 kWe electrical loss, 5 kWe parasitic", "fsp2010icd",
         FSPInputs(delivered_power_kwe=40, pmad_efficiency=40 / 43, fixed_parasitic_kwe=5.0,
                   efficiency_mode="fixed", conversion_efficiency=0.26,
                   hot_temperature_k=778, cold_temperature_k=425),
         {"converter_output_kwe": 48.0, "reactor_thermal_kwt": 186.0,
          "converter_heat_rejection_kwt": 140.0}, {"converter_output_kwe"}),
        ("S2", "2025 HP-Stirling: PMAD 83 %", "mason2025trades",
         FSPInputs(delivered_power_kwe=40, pmad_efficiency=0.83, fixed_parasitic_kwe=p_fix_2025,
                   efficiency_mode="fixed", conversion_efficiency=50.4 / 175,
                   hot_temperature_k=1050, cold_temperature_k=420, sink_temperature_k=270),
         {"converter_output_kwe": 50.4, "converter_heat_rejection_kwt": 119.0}, set()),
        ("S3", "2025 GC-Brayton: PMAD 90 %", "mason2025trades",
         FSPInputs(delivered_power_kwe=40, pmad_efficiency=0.90, fixed_parasitic_kwe=p_fix_2025,
                   efficiency_mode="fixed", conversion_efficiency=46.7 / 230,
                   hot_temperature_k=1100, cold_temperature_k=390, sink_temperature_k=270),
         {"converter_output_kwe": 46.7, "converter_heat_rejection_kwt": 180.0}, set()),
    ]


# Scaling: Table V of mason2025trades (reactor kWt and radiator m^2 vs net power).
SCALING = {
    "HP-Stirling": {"pmad": 0.83, "eta": 50.4 / 175,
                    "kwt": {10: 49, 20: 90, 40: 175, 80: 343, 100: 429},
                    "m2": {10: 40, 20: 80, 40: 140, 80: 279, 100: 351}},
    "GC-Brayton": {"pmad": 0.90, "eta": 46.7 / 230,
                   "kwt": {10: 72, 20: 120, 40: 230, 80: 421, 100: 517},
                   "m2": {10: 57, 20: 93, 40: 180, 80: 317, 100: 380}},
}


def scaling_results():
    """Predict Table V from the 40 kWe point (same eta and PMAD efficiency),
    with no fixed load and with a 2 kWe fixed load (an interpretation)."""
    out = {}
    for name, d in SCALING.items():
        rows = []
        preds = {}
        for pfix in (0.0, 2.0):
            base = FSPInputs(pmad_efficiency=d["pmad"], fixed_parasitic_kwe=pfix,
                             efficiency_mode="fixed", conversion_efficiency=d["eta"],
                             hot_temperature_k=1100, cold_temperature_k=420)
            r40 = size_fsp(replace(base, delivered_power_kwe=40))
            for p in d["kwt"]:
                r = size_fsp(replace(base, delivered_power_kwe=p))
                # radiator scales with heat load; anchored at the published 40 kWe area
                preds[(pfix, p)] = (r.reactor_thermal_kwt,
                                    d["m2"][40] * r.converter_heat_rejection_kwt
                                    / r40.converter_heat_rejection_kwt)
        for p, kwt_pub in d["kwt"].items():
            rows.append({"kwe": p, "kwt_pub": kwt_pub,
                         "kwt_pred0": preds[(0.0, p)][0], "kwt_pred2": preds[(2.0, p)][0],
                         "m2_pub": d["m2"][p], "m2_pred2": preds[(2.0, p)][1],
                         # efficiency NASA's numbers imply at this size (2 kWe fixed load)
                         "eta_implied": (p / d["pmad"] + 2.0) / kwt_pub})
        out[name] = rows
    return out


# Implied fraction of Carnot: f = eta / (1 - Tc/Th) for each published point.
CARNOT_POINTS = [
    ("2010 reference, effective cycle temps", "fsp2010icd", 0.26, 778, 425),
    ("2025 HP-Stirling, hot/cold-end temps", "mason2025trades", 50.4 / 175, 1050, 420),
    ("2025 GC-Brayton, turbine inlet/cold end", "mason2025trades", 46.7 / 230, 1100, 390),
    ("Deployable concept, Stirling (700/150 C)", "oleson2022deployable", 0.261, 973.15, 423.15),
    ("TDU test, max power (NaK 875 K, water 350 K)", "briggs2016tdu", 0.233, 875, 350),
    ("SSE CTPC test, 12.5 kWe (1050 K / 525 K)", "dhar1999sse", 0.22, 1050, 525),
]


def pct(pred, pub):
    return (pred - pub) / pub


def run():
    results = {"radiator": [], "system": [], "scaling": {}, "carnot": []}
    b1 = RADIATOR_CASES[0]
    eta_cal = calibrate_eta_p(b1)
    results["eta_p_calibrated"] = eta_cal

    for c in RADIATOR_CASES:
        ideal = radiator_area(c, 1.0)
        pred = radiator_area(c, eta_cal)
        role = "calibration" if c.key == "B1" else "prediction"
        results["radiator"].append({
            "key": c.key, "label": c.label, "source": c.source, "role": role,
            "q_kwt": c.q_kwt, "t_in": c.t_in, "t_out": c.t_out, "t_sink": c.t_sink,
            "area_pub": c.area_pub_m2, "area_ideal": ideal, "area_pred": pred,
            "err_ideal": pct(ideal, c.area_pub_m2), "err_pred": pct(pred, c.area_pub_m2),
            "eta_p_implied": calibrate_eta_p(c),
            "pass": role == "calibration" or abs(pct(pred, c.area_pub_m2)) <= TOL_RADIATOR,
            "assumed": c.assumed})

    for pfix in (2.0, 0.0):
        for key, label, src, inp, pub, ident in system_cases(pfix):
            if key == "S1" and pfix == 0.0:
                continue
            r = size_fsp(inp)
            for metric, val in pub.items():
                pred = getattr(r, metric)
                is_id = metric in ident
                results["system"].append({
                    "key": key, "label": label, "source": src, "metric": metric,
                    "pfix": inp.fixed_parasitic_kwe if key != "S1" else None,
                    "pub": val, "pred": pred, "err": pct(pred, val), "identity": is_id,
                    "pass": True if is_id else abs(pct(pred, val)) <= TOL_POWER})

    # Electronics radiator, information only (no acceptance criterion).
    # Published: controller radiator + user-interface radiator [mason2025trades].
    results["electronics"] = []
    for (key, label, src, inp, pub, ident), published in zip(system_cases(2.0)[1:], (18 + 7, 6 + 7)):
        r = size_fsp(inp)
        results["electronics"].append({"key": key, "pub": published,
                                       "pred": r.electronics_radiator_area_m2,
                                       "q_kw": r.electronics_heat_kwt, "t_sink": inp.sink_temperature_k})

    # Sensitivity of the B2/B3 predictions to the (unstated) sizing sink.
    results["sink_sensitivity"] = {}
    for c in RADIATOR_CASES[1:3]:
        c250 = replace(c, t_sink=250.0)
        results["sink_sensitivity"][c.key] = pct(radiator_area(c250, eta_cal), c.area_pub_m2)

    results["scaling"] = scaling_results()
    for label, src, eta, th, tc in CARNOT_POINTS:
        results["carnot"].append({"label": label, "source": src, "eta": eta, "th": th, "tc": tc,
                                  "carnot": 1 - tc / th, "f": eta / (1 - tc / th)})
    return results


def write_latex(res, path="docs/benchmark_tables.tex"):
    """Generate the LaTeX tables the theory document \\input{}s, so the
    document always shows what the code actually computes."""
    rad = {r["key"]: r for r in res["radiator"]}
    el = {r["key"]: r for r in res["electronics"]}
    mac = {
        "etapcal": f'{res["eta_p_calibrated"]:.3f}',
        "errBtwo": f'{100*rad["B2"]["err_pred"]:+.1f}',
        "errBthree": f'{100*rad["B3"]["err_pred"]:+.1f}',
        "errBfour": f'{100*rad["B4"]["err_pred"]:+.1f}',
        "idealmin": f'{-100*max(r["err_ideal"] for r in res["radiator"]):.0f}',
        "idealmax": f'{-100*min(r["err_ideal"] for r in res["radiator"]):.0f}',
        "sinkBtwo": f'{100*res["sink_sensitivity"]["B2"]:+.1f}',
        "sinkBthree": f'{100*res["sink_sensitivity"]["B3"]:+.1f}',
        "elecStirling": f'{el["S2"]["pred"]:.0f}', "elecBrayton": f'{el["S3"]["pred"]:.0f}',
        "elecStirlingPub": f'{el["S2"]["pub"]:.0f}', "elecBraytonPub": f'{el["S3"]["pub"]:.0f}',
    }
    header = ["% AUTO-GENERATED by benchmarks.py -- do not edit by hand."]
    header += [r"\newcommand{\%s}{%s}" % (k, v) for k, v in mac.items()]
    blocks = {}
    L = []
    L.append(r"\begin{table}[ht]\centering\small")
    L.append(r"\caption{Main-radiator benchmark. Areas are total radiating area (both faces) "
             r"including margin. B1 is used to calibrate $\eta_p$; B2--B4, not used in the "
             r"calibration, are predicted with $\eta_p=\etapcal$ held fixed.}\label{tab:radbench}")
    L.append(r"\resizebox{\textwidth}{!}{\begin{tabular}{llrrrrrrr}\toprule")
    L.append(r"Case & Role & $Q$ (kW) & $T_\text{in}/T_\text{out}$ (K) & $T_s$ (K) & "
             r"Published (m$^2$) & Ideal $\eta_p{=}1$ & Calibrated & Implied $\eta_p$\\\midrule")
    for r in res["radiator"]:
        L.append(rf"{r['key']}~\cite{{{r['source'].split(';')[0].strip()}}} & {r['role']} & "
                 rf"{r['q_kwt']:.1f} & {r['t_in']:.0f}/{r['t_out']:.0f} & {r['t_sink']:.0f} & "
                 rf"{r['area_pub']:.1f} & {r['area_ideal']:.1f} ({100*r['err_ideal']:+.0f}\%) & "
                 rf"{r['area_pred']:.1f} ({100*r['err_pred']:+.1f}\%) & {r['eta_p_implied']:.3f}\\")
    L.append(r"\bottomrule\end{tabular}}\end{table}")
    blocks["tableradbench"], L = L, []

    L.append(r"\begin{table}[ht]\centering\small")
    L.append(r"\caption{Electrical chain and first-law benchmark (acceptance $\pm5\%$). "
             r"$P_\text{fix}$ is the fixed load assumed for the 2025 designs (an interpretation, "
             r"see text). Rows marked $\dagger$ are identities by construction and are not counted "
             r"as agreement.}\label{tab:sysbench}")
    L.append(r"\begin{tabular}{llcrrr}\toprule Case & Quantity & $P_\text{fix}$ (kWe) & Published & TLECUIL & Difference\\\midrule")
    names = {"converter_output_kwe": "Gross converter output (kWe)",
             "reactor_thermal_kwt": "Reactor thermal power (kWt)",
             "converter_heat_rejection_kwt": "Main radiator heat load (kWt)"}
    for r in res["system"]:
        pf = "--" if r["pfix"] is None else f"{r['pfix']:.0f}"
        dag = r"$^\dagger$" if r["identity"] else ""
        L.append(rf"{r['key']}~\cite{{{r['source']}}} & {names[r['metric']]}{dag} & {pf} & {r['pub']:.1f} & "
                 rf"{r['pred']:.1f} & {100*r['err']:+.1f}\%\\")
    L.append(r"\bottomrule\end{tabular}\end{table}")
    blocks["tablesysbench"], L = L, []

    L.append(r"\begin{table}[ht]\centering\small")
    L.append(r"\caption{Scaling benchmark against Table~V of \cite{mason2025trades}. TLECUIL is "
             r"anchored only at 40\,kWe (same $\eta$ and PMAD efficiency) and predicts the other "
             r"power levels with no fixed load and with a 2\,kWe fixed load. Implied $\eta$ is the "
             r"efficiency NASA's numbers imply at each size.}\label{tab:scalebench}")
    L.append(r"\resizebox{\textwidth}{!}{\begin{tabular}{llrrrrrrrrr}\toprule Concept & kWe & "
             r"kWt pub. & kWt ($P_\text{fix}{=}0$) & Diff. & kWt ($P_\text{fix}{=}2$) & Diff. & "
             r"Implied $\eta$ & m$^2$ pub. & m$^2$ ($P_\text{fix}{=}2$) & Diff.\\\midrule")
    for name, rows in res["scaling"].items():
        for i, r in enumerate(rows):
            L.append(rf"{name if i == 0 else ''} & {r['kwe']} & {r['kwt_pub']} & {r['kwt_pred0']:.0f} & "
                     rf"{100*pct(r['kwt_pred0'], r['kwt_pub']):+.1f}\% & {r['kwt_pred2']:.0f} & "
                     rf"{100*pct(r['kwt_pred2'], r['kwt_pub']):+.1f}\% & {r['eta_implied']:.3f} & "
                     rf"{r['m2_pub']} & {r['m2_pred2']:.0f} & {100*pct(r['m2_pred2'], r['m2_pub']):+.1f}\%\\")
        L.append(r"\addlinespace")
    L.append(r"\bottomrule\end{tabular}}\end{table}")
    blocks["tablescalebench"], L = L, []

    L.append(r"\begin{table}[ht]\centering\small")
    L.append(r"\caption{Fraction of Carnot implied by published efficiencies, "
             r"$f_C=\eta/(1-T_c/T_h)$. The value depends on how $T_h$ and $T_c$ are "
             r"defined.}\label{tab:carnot}")
    L.append(r"\begin{tabular}{lrrrrr}\toprule Source & $\eta$ & $T_h$ (K) & $T_c$ (K) & $\eta_C$ & $f_C$\\\midrule")
    for r in res["carnot"]:
        L.append(rf"{r['label']}~\cite{{{r['source']}}} & {r['eta']:.3f} & {r['th']} & {r['tc']} & "
                 rf"{r['carnot']:.3f} & {r['f']:.2f}\\")
    L.append(r"\bottomrule\end{tabular}\end{table}")
    blocks["tablecarnot"] = L
    out = header[:]
    for name, lines in blocks.items():
        out.append(r"\newcommand{\%s}{%%" % name)
        out.extend(lines)
        out.append("}")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return path


def main():
    res = run()
    print(f"Panel effectiveness calibrated on B1 only: eta_p = {res['eta_p_calibrated']:.3f}\n")
    print("RADIATOR (total radiating area incl. margin, m^2)")
    print(f"{'case':4} {'role':12} {'pub':>7} {'ideal':>14} {'calibrated':>16} {'eta_p impl.':>11}  pass")
    for r in res["radiator"]:
        print(f"{r['key']:4} {r['role']:12} {r['area_pub']:7.1f} "
              f"{r['area_ideal']:7.1f} ({100*r['err_ideal']:+4.0f}%) "
              f"{r['area_pred']:8.1f} ({100*r['err_pred']:+5.1f}%) {r['eta_p_implied']:11.3f}  {r['pass']}")
    print("\nSYSTEM (identity rows are not agreement)")
    for r in res["system"]:
        tag = "identity" if r["identity"] else r["pass"]
        print(f"{r['key']:4} Pfix={r['pfix']} {r['metric']:30} pub {r['pub']:7.1f}  pred {r['pred']:7.1f}  "
              f"{100*r['err']:+5.1f}%  {tag}")
    print("\nELECTRONICS RADIATOR (info only)")
    for r in res["electronics"]:
        print(f"{r['key']:4} pub {r['pub']} m2  pred {r['pred']:.1f} m2  (Q {r['q_kw']:.1f} kW, sink {r['t_sink']} K)")
    print("\nSINK SENSITIVITY (B2/B3 at 250 K instead of 270 K):", {k: f"{100*v:+.1f}%" for k, v in res["sink_sensitivity"].items()})
    print("\nSCALING (reactor kWt; radiator m^2)")
    for name, rows in res["scaling"].items():
        for r in rows:
            print(f"{name:12} {r['kwe']:4} kWe  kWt {r['kwt_pub']:4} vs {r['kwt_pred0']:6.1f} "
                  f"({100*pct(r['kwt_pred0'], r['kwt_pub']):+5.1f}%) / Pfix=2 {r['kwt_pred2']:6.1f} "
                  f"({100*pct(r['kwt_pred2'], r['kwt_pub']):+5.1f}%)   m2 {r['m2_pub']:4} vs "
                  f"{r['m2_pred2']:6.1f} ({100*pct(r['m2_pred2'], r['m2_pub']):+5.1f}%)  eta_impl {r['eta_implied']:.3f}")
    print("\nIMPLIED FRACTION OF CARNOT")
    for r in res["carnot"]:
        print(f"  {r['label']:48} f = {r['f']:.2f}")
    with open("benchmark_results.json", "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print("\nwrote benchmark_results.json and", write_latex(res))


if __name__ == "__main__":
    main()
