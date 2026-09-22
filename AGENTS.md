# Repository instructions

## Manuscript and integrity

- The canonical source is `deterministic_von_neumann_fair_extractor.tex`; its PDF is tracked.
- Preserve mathematical statements, hypotheses, proof content, authorship, notation, theorem numbering, and stable labels unless explicitly instructed to change them.
- The manuscript has one author: Christopher D. Long. Do not add coauthors without explicit instruction.
- Flag suspected mathematical errors; do not silently repair, weaken, or broaden claims.
- Keep mathematical changes separate from migration, formatting, and repository-maintenance changes.
- Compile with `latexmk -pdf -interaction=nonstopmode -halt-on-error deterministic_von_neumann_fair_extractor.tex` and check undefined citations/references and layout warnings before committing a new PDF.
- Do not commit LaTeX auxiliary files or build caches.

## Exact verification

- Exact verification code, checked outputs, and finite certificates belong in `scripts/`.
- The canonical verifier is `scripts/verify_glazer_critical.py`. Its SHA-256 digest recorded by the manuscript is `07801d62aedf91b87e122bd96f53a9bd450d64dd4e265caaff364f3da0dced61`; do not change the verifier without deliberately updating and re-auditing that manuscript statement.
- Required witnesses are in `scripts/certificates/finite_blocks.json`.
- Run the verifier both constructively and in independent certificate-checking mode, and run `scripts/supplementary_checks.py`.
- Do not silently regenerate recorded outputs merely to make checks pass.
- The finite computations certify the stated finite blocks and implementation identities; the all-dyadic theorem is the mathematical argument in the manuscript, not an extrapolation from testing.

## Formalization

- No Lean project is currently part of this repository.
- `FORMALIZATION_PLAN.md` and `FORMALIZATION_LEDGER.csv` are planning documents, not a formal verification claim.
- If formalization is explicitly started, follow the current `long-mathematics` conventions: root Lake project, named module directory, root umbrella module, explicit status ledger, and dependency/axiom audits.
- Never claim formal verification before statement correspondence and transitive dependency audits are complete.

## Git workflow

Use a feature branch, pull request, successful applicable CI, and squash merge. Never force-push `main` or overwrite unique source material.

## Documentation

Keep the repository a current paper companion rather than a development archive. The README should describe the theorem, manuscript links, exact checks, scope limitations, formalization status, citation, and license. Do not add obsolete drafts, local machine paths, or redundant build logs.
