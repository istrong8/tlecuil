## Change
<!-- What changed and why -->

## Checklist (V&V plan Secs. 6-8)
- [ ] Tests pass (`python -m pytest -q`)
- [ ] If physics changed: Theory document updated, PDF rebuilt, `python docs/sync_refs.py --fix` run
- [ ] If results changed: version bumped and CHANGELOG entry added
- [ ] Same change made in `web/tlecuil.html`; parity check passes
- [ ] Benchmark tables regenerated (`python benchmarks.py`)
- [ ] Independent reviewer assigned if equations changed
