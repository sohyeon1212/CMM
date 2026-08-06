# CLAUDE.md

@AGENTS.md

## Claude Code specifics

- The scenario pipelines in `docs/scenarios/` and the function reference in
  `docs/agent-reference.md` are loaded on demand. Read the one file you need rather than
  globbing the directory.
- Long solves (FVA on genome-scale models, OptKnock/RobustKnock, large sampling runs) can run
  for minutes. Run them with a raised `timeout` or in the background rather than killing and
  retrying with smaller parameters — a silently shrunk `n_steps` or sample count changes the
  science.
- Driving the Qt app headlessly (only needed when testing the GUI itself):
  ```bash
  QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras python -m cmm.app
  QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest tests/test_app_smoke.py -q
  ```
- Use `uv run --frozen --all-extras` for anything that must match the locked publication
  environment; the bare `.venv/bin/python` does not have the test and lint tooling installed.
