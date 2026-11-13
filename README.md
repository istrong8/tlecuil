# TLECUIL: A First-Order Fission Surface Power Sizer

*Tlecuilli* (Nahuatl) means the three-stone hearth. The three stones stand for the
three subsystems this code balances: the reactor, power conversion and heat rejection.

TLECUIL is an open, independent, first-order steady-state thermal-electrical
balance for fission surface power (FSP) systems. It has no power limit. Every equation
is derived and cited in `docs/tlecuil_theory.pdf`. The code is verified by 16 tests
and benchmarked against four published NASA designs.

| Path | What it is |
|---|---|
| `tlecuil.py` | The engine (reference implementation). Comments cite the theory equations and the sources. |
| `benchmarks.py` | Comparison with NASA designs (2010 reference, 2025 GRD Stirling and Brayton, deployable 40 kWe). Regenerates `docs/benchmark_tables.tex`. |
| `test_tlecuil.py` | Verification and benchmark tests. |
| `web/tlecuil.html` | Interactive page. Same engine in JavaScript, with presets for the NASA designs. |
| `check_parity.js` | Proves that the HTML and Python engines give identical numbers. |
| `docs/tlecuil_theory.tex` | Theory for undergraduates: primer, derivations, worked example, exercises, benchmarks, EZ FSP Sizer comparison, Moon vs Mars. |
| `docs/tlecuil_vv_plan.tex` | Software quality assurance and V&V plan (NQA-1-2024 Subpart 2.7 structure). |
| `docs/make_figures.py` | Regenerates the figures from the code. |
| `docs/sync_refs.py` | Keeps equation numbers in code comments in step with the PDF. |
| `docs/references.bib` | Bibliography (biblatex + biber). |
| `.vscode/settings.json` | LaTeX Workshop recipe for VS Code. |

## Quick start

```bash
pip install numpy scipy matplotlib pytest
python tlecuil.py                    # generic 40 kWe case
python tlecuil.py --power 100        # any power
python benchmarks.py                 # NASA comparison; rewrites docs/benchmark_tables.tex
python -m pytest -q                  # 16 tests
```

## Building the PDF in VS Code

1. Install a TeX distribution: TeX Live (full) or MiKTeX. It must include `biblatex`,
   `biblatex-ieee` and `biber`.
2. Install the **LaTeX Workshop** extension.
3. Open `docs/tlecuil_theory.tex` and press **Ctrl+Alt+B** (Cmd+Option+B on Mac). The
   recipe in `.vscode/settings.json` runs `latexmk` with biber.
   - Command line alternative: `cd docs && latexmk -pdf tlecuil_theory.tex`
4. If you change the code, run `python benchmarks.py` before rebuilding, so the
   benchmark tables match the code.

## Keeping Python, HTML and PDF in sync

1. Change the equation in `tlecuil.py` and add a test.
2. Make the same change in the `ENGINE-START … ENGINE-END` block of `web/tlecuil.html`.
3. Run these commands. The parity check must report a difference of 0.

   ```bash
   python tlecuil.py --export-json parity_cases.json
   node check_parity.js web/tlecuil.html parity_cases.json
   ```
4. After rebuilding the PDF, run `python docs/sync_refs.py --fix`. It updates every `eq. X.Y [label]` reference in the code to match the PDF, and `--check` is run in CI.

## What it does not do

TLECUIL does not model mass, shielding, dose, reactor physics, transients or
size-dependent converter efficiency. The benchmarks are model-to-model comparisons
with published design studies, not hardware validation. See Theory Sec. 14.

## Quality assurance

TLECUIL is developed under a graded software QA and V&V plan
(`docs/tlecuil_vv_plan.pdf`), structured after ASME NQA-1-2024 Part II Subpart 2.7, so it
can later be qualified by a formal QA program. CI (`.github/workflows/ci.yml`) runs on
every change: the tests, a check that the benchmark tables are current, the
JavaScript/Python parity check, the PDF build, and the reference sync check. Report
errors with the **Problem report** issue template.

## License

BSD 3-Clause. See `LICENSE`. To cite TLECUIL, use `CITATION.cff`.
