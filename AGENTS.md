# Repository Guidelines

This repository mixes product code, technical documentation, literature intake, and report deliverables. Keep those concerns separated.

## Compatibility Constraints

- Keep `sx_dashboard.py` at the repository root. The current web app entrypoint depends on that path.
- Keep `sx_simulator_app.py` at the repository root unless every referenced command and document is updated together.
- Keep `sx_simulator/` as the runtime package root.
- Keep these markdown paths stable because the web app loads them directly:
  - `docs/user_manual.md`
  - `docs/glossary.md`
  - `docs/references.md`
  - `docs/validation_history.md`
- Keep `docs/literature/` as the literature intake root. Recent local literature workflow and handoff documents depend on that path.

## Placement Rules

- New runtime or library code goes in `sx_simulator/`.
- New app entry scripts should prefer `scripts/` unless they must remain root-level for compatibility.
- New one-off analysis or comparison scripts go in `scripts/analysis/`.
- New verification helpers or report-style validation scripts go in `scripts/verification/` unless they intentionally replace a root compatibility script.
- User-facing or architecture documentation goes in `docs/`.
- Literature PDFs and literature handoff notes go in `docs/literature/`.
- External-facing report decks, briefing packages, and draft workspaces go in `deliverables/`.
- Temporary renders, extracted text, caches, and local scratch outputs go in `tmp/` and should stay out of git.

## Naming Rules

- Prefer English snake_case for code files and directories.
- Preserve Korean filenames for business-facing documents when that matches existing stakeholder usage.
- Prefer directory-based grouping over long filenames that encode category, audience, and lifecycle all at once.

## Before Moving Anything

- Search for path references with `rg`.
- Preserve or replace compatibility paths used by the web app, scripts, or literature workflow.
- If a file is cited by many docs with line-specific links, prefer leaving it in place and document the rule instead of moving it.
