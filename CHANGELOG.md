# Changelog

All notable changes are recorded here. A change that alters any numerical
result increments at least the MINOR version (see docs/tlecuil_vv_plan.pdf, Sec. 8).

## [0.2.1] - 2026-11-13
### Changed
- Review IR-001 relabelled as AI-assisted (performed by a separate AI model
  instance), not independent human review. IR-002, an independent human review,
  added to the V&V plan as planned.
- Author name updated to Irais Strong.
- Interactive page published at https://istrong8.github.io/tlecuil/.
No numerical results changed.

## [0.2.0] - 2026-11-11
### Fixed
- PR-001: v0.1 divided by the number of radiating sides twice. The "emitting area"
  it reported was half the true radiating area for two-sided panels
  (non-conservative). Regression test added.
- PR-002: the closed-form radiator integral lost precision for sink temperatures
  near 0 K. Series branch added; tested against quadrature from 1e-9 K to 290 K.
### Added
- Renamed to TLECUIL. Fixed parasitic loads, reactor heat losses, panel
  effectiveness, area margin.
- Benchmarks against four published NASA designs (benchmarks.py).
- Theory document rewritten for undergraduates: primer, figures, worked example,
  exercises, glossary; biblatex references.
- Software QA and V&V plan (docs/tlecuil_vv_plan.tex).
- docs/sync_refs.py: keeps code comments in step with equation numbers.
- Interactive page with NASA design presets; JavaScript/Python parity check.
- AI-assisted technical review IR-001; findings dispositioned.

## [0.1.0] - 2026-11-11
- First draft (fsp_sizer.py). Superseded; see PR-001.
