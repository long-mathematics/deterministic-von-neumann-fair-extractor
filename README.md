# A Critical Linear Bound with Logarithmic Savings for Glazer's Time-Limited Fair-Coin Problem

## Abstract

Let a coin have unknown bias `p in (0,1)`, and let `L` be the length of its maximal initial constant run. This repository accompanies a deterministic, bias-independent fair-bit extractor satisfying

```text
T(L) <= ceil(C_* L)
T(L) <= C_* L - log_5(L) + 6
C_* = 1 + 2 log(phi)/log(5) = 1.597987435665...
phi = (1 + sqrt(5))/2.
```

The construction separately balances every dyadic block at common horizon `2N`. For `N <= L < 2N`, with `N >= 2` a power of two, the exact blockwise guarantee is

```text
s_N = max(0, floor(log_5(N)) - 3)
T(L) <= min(2N, ceil(C_* L) - s_N).
```

The argument is uniform for dyadic `N >= 128`; the six smaller blocks have exact integer certificates. The leading coefficient `C_*` is optimal for all separately balanced horizon-`2N` rules, without assuming complement symmetry, and the lower bound persists with `o(N)` additional horizon. The manuscript also bounds the possible logarithmic saving, extends the fairness statement to exchangeable binary sources with almost surely finite initial run, and characterizes why powers of two are the admissible block lengths at horizon `2N`.

## Manuscript and source

- [Research draft PDF](deterministic_von_neumann_fair_extractor.pdf)
- [LaTeX source](deterministic_von_neumann_fair_extractor.tex)
- [Exact verifier](scripts/verify_glazer_critical.py)
- [Finite packet certificates](scripts/certificates/finite_blocks.json)
- [Recorded verification output](scripts/verification_output.txt)
- [Supplementary exact checks](scripts/supplementary_checks.py)
- [Formalization plan](FORMALIZATION_PLAN.md)
- [Formalization ledger](FORMALIZATION_LEDGER.csv)

The manuscript is a research draft dated September 20, 2026. It has been model-assisted and internally audited, but is not yet independently peer reviewed or formally verified.

## What is improved over the earlier critical-bound draft

The leading constant remains

```text
C_* = 1 + 2 log(phi)/log(5).
```

The revised theorem strengthens the earlier ceiling bound in several directions:

- a logarithmic saving `T(L) <= C_* L - log_5(L) + 6` while preserving separate fairness within each dyadic block;
- optimality of the leading coefficient for the full separately balanced horizon-`2N` class, not only complement-symmetric packet rules;
- persistence of that lower bound under `o(N)` additional horizon;
- an exact deadline inequality restricting second-order logarithmic improvements;
- extension from fixed-bias i.i.d. inputs to exchangeable binary sources with almost surely finite initial run;
- a necessity-and-sufficiency explanation for dyadic block lengths at horizon `2N`, together with a two-output-bit obstruction.

The paper does **not** claim unrestricted optimality among all globally fair extractors when imbalance may cancel between different dyadic blocks.

## Reproducible exact checks

Python 3.10+ is sufficient; the verification code uses only the standard library. Do not use Python's `-O` flag because assertions are part of the checker.

```sh
python3 scripts/verify_glazer_critical.py
python3 scripts/verify_glazer_critical.py \
  --check-certificate scripts/certificates/finite_blocks.json
python3 scripts/verify_glazer_critical.py --extra 128 256 512 1024 \
  --export /tmp/finite_blocks.json
python3 scripts/supplementary_checks.py
```

The verifier's SHA-256 digest is

```text
07801d62aedf91b87e122bd96f53a9bd450d64dd4e265caaff364f3da0dced61
```

which is the digest recorded in the manuscript. The extended run checks 471,442 packet sites for `N=1024`; all supplied checks use exact integer or rational arithmetic.

Compile the manuscript with:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error deterministic_von_neumann_fair_extractor.tex
```

## Formalization status

No Lean formalization is claimed. `FORMALIZATION_PLAN.md` records a proposed decomposition of the proof, and `FORMALIZATION_LEDGER.csv` inventories the numbered mathematical statements and auxiliary obligations. They are planning documents only.

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). The draft does not yet have an arXiv identifier or DOI.

## License

Code and repository infrastructure are provided under the [MIT License](LICENSE). The mathematical manuscript remains attributable to Christopher D. Long as stated in the paper.
