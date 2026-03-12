# Repository Organization

## Goal

Keep product code, research material, validation work, and stakeholder deliverables separate so the repository stays navigable without breaking the current web app or local literature workflow.

## Stable Paths

These paths are intentionally stable and should not be moved casually:

- `sx_dashboard.py`
- `sx_simulator_app.py`
- `sx_simulator/`
- `docs/user_manual.md`
- `docs/glossary.md`
- `docs/references.md`
- `docs/validation_history.md`
- `docs/literature/`

## Folder Roles

- Root: only entrypoints, top-level project docs, and essential config.
- `sx_simulator/`: runtime engine and app-facing service code.
- `tests/`: automated tests that are expected to fail on regression.
- `scripts/`: manual analysis or helper scripts.
- `docs/`: documentation and review records.
- `docs/analysis/`: technical notes or resolved summaries that are not part of the web app.
- `docs/architecture/`: repository and code-structure guidance.
- `docs/literature/`: literature intake, categorized papers, and handoff notes.
- `deliverables/`: presentations, report packages, and stakeholder-facing workspaces.
- `tmp/`: scratch outputs, extracted text, rendered images, and local caches.

## Current Decisions

- The web app still reads selected docs from `docs/` directly, so those files remain flat for compatibility.
- `validation_test.py` and `test_verification.py` remain at the root because many historical documents refer to those exact paths.
- Literature stays under `docs/literature/` because the current local literature workflow and handoff notes already use that structure.

## Ongoing Rule

When adding a new file, decide first whether it is:

1. runtime code,
2. automated verification,
3. manual analysis,
4. documentation,
5. literature,
6. stakeholder deliverable, or
7. temporary output.

Then place it in the matching directory instead of the repository root.
