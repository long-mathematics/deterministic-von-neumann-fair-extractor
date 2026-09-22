# Formalization plan

## Status and intended scope

This is a plan, not a Lean development or a completion report. No Lean theorem has been checked in this package. The mathematical source of truth is `deterministic_von_neumann_fair_extractor.tex`; the finite certificate data are in `scripts/certificates/finite_blocks.json`.

The intended scope is the full paper: the logarithmically improved extractor, all auxiliary construction lemmas, the unrestricted-signing block obstruction, the second-order bound, the sublinear-extra-horizon extension, exchangeable-source fairness, and the dyadic and output-cardinality results. The original ceiling bound alone is not full coverage.

There is no apparent need for a deep external theorem outside the elementary material proved in the manuscript. This is a favorable formalization target, but the custom rounding argument and correspondence between finite signings and an infinite stopping rule require substantive work. Feasibility has not been established by a trial Lean build.

## Statement interfaces to settle first

Use an infinite input `Nat -> Bool` only at the outermost layer. Initially formalize finite words, a fixed zero-count statistic, run classes, and prefix-determined sign functions. The principal counting theorem should state equality of signed counts for every zero count at a fixed horizon.

Define a stopping rule with prefix measurability, deterministic output, and a bound on every continuation with finite initial run length. The bound must not mean an expectation, a high-probability deadline, or merely termination almost surely under one bias. Handle `L=1` explicitly and leave the two constant streams outside the pathwise finite-termination domain.

Keep the three fairness notions separate:

1. Unconditional fairness for every fixed bias in `(0,1)`.
2. Equal output masses separately within a dyadic block.
3. Equal output masses conditional on a dyadic block and the padded zero count.

The manuscript does not claim fairness conditional on each individual run length or independence from the actual stopping time.

For the lower bounds, do not assume complement symmetry. The crucial hypotheses are separate block balancing and the horizon restriction. For the robust version, the horizon is `2N+j_N`, with nonnegative integer `j_N` and `j_N/N -> 0` along dyadic `N`.

## Suggested dependency order

### 1. Finite packet and type-count layer

Use `Nat.choose`, finite words, integer signed counts, and univariate polynomials. Prove packet realization by selecting exactly `(capacity + coefficient)/2` words in each weight class. Formalize both directions of type balance. Prove that a per-level deadline makes the sign constant on unread suffix cylinders, even for an adaptive stopping algorithm padded to that deadline.

### 2. Constants and profiles

Define `phi`, `beta`, `alpha`, `C_*`, and the integer shift. Prove positivity, irrationality, the critical power identities, and the elementary rational bounds. Treat integer/natural casts and subtraction bounds explicitly rather than relying on truncated subtraction silently.

For the six finite blocks, certify `floor(beta*m)` algebraically. The equivalence is

```text
floor(beta*m) = r  <=>  5^r < phi^(2m) < 5^(r+1),
```

for positive integer `m`. Express `2*phi^(2m)` as `A+B*sqrt(5)` by the proved Fibonacci recurrence, and reduce ordered comparisons to integer-square inequalities. This avoids making a verified transcendental-series evaluator a prerequisite for the finite witnesses. The decimal interval remains an optional independently verified numerical corollary.

### 3. Folding and finite binomial estimates

Represent polynomials by finite coefficient data before lifting them to an abstract polynomial interface. Prove the one-sign recursion and telescoping identity, then the exact folded-mass formula with the stated monotonicity hypotheses.

The new atom estimate can be proved entirely as a finite mixture identity plus weighted Cauchy–Schwarz. The level-zero tail estimate uses the finite second moment `n/4`. No Fourier inversion, Gaussian integral, Hoeffding theorem, de Finetti theorem, or infinite negative-binomial series is needed.

Keep these finite-sum lemmas independent of general probability measures. A later adapter can connect them to Mathlib's binomial distribution if useful. The current binomial API exposes explicit singleton masses and finite-sum integral formulas, but exact APIs should be inspected again after the toolchain is pinned.

### 4. Integer lattice rounding

This is the highest-risk custom proof-engineering component. Use an explicit finite dependent type of packet sites `(i,r)` and ensure every child index is in range. Prove the unitriangular boundary basis and the parent/child kernel basis. Construct the parity-correct integral reference point.

The rounding lemma must retain all structural hypotheses stated in the revision: suffix drops in `{0,1,2}`, the final zero-suffix segment, the absence of early overlap of the two boundary chains, and vanishing of the center on the last two levels. Formalize the lower-boundary prescription, upper-chain induction, and the final shared variable separately. In particular, prove rather than assume that the preceding upper parameter is an even integer of absolute value at most one and hence zero.

Then prove the five-unit interior error and the stronger two-unit zero-tail error. These are distinct bounds and are both used in the all-block theorem.

### 5. Finite certificates in the proof assistant

Separate generation from checking. The external Python generator may produce candidate arrays, but Lean should check the arrays through a proved checker or explicit proof-producing arithmetic. A successful Python run is not an axiom establishing their validity.

Only the six blocks below 128 are required. Formalize the algebraic floor witnesses and coefficient-by-coefficient identities for those blocks. Large exploratory blocks need not be replayed in the kernel to obtain the theorem. Prefer kernel-checkable computation or proof-producing reflection; any extension to the trusted computation base must be documented, not silently substituted for ordinary kernel checking.

### 6. Uniform construction and infinite-source wrapper

Combine the norm, reserve, and rounding results for every dyadic `N>=128`, and the six checked witnesses below the threshold. Choose the per-block signings deterministically; the paper's lexicographic realization supplies an explicit option. Prove the pathwise deadline and both advertised inequalities before introducing an input measure.

Finally, prove cylinder measurability and disjointness of dyadic run events. For product Bernoulli inputs, the two constant streams have measure zero. Countable additivity transfers blockwise balance to global fairness. For exchangeable sources, use only equality of finite-word probabilities within each zero-count class and the explicit assumption that the initial run is finite almost surely.

### 7. Lower bounds and extensions

The horizon-`2N` obstruction is a finite polynomial theorem. Separate zero-run and one-run contributions by total zero count, prove the zero-run polynomial is the constant `+1` or `-1`, and derive the exact deadline inequality. The leading and logarithmic obstructions then use elementary exponential/logarithm estimates.

The `2N+o(N)` extension adds a degree bound, finite Lagrange interpolation on a fixed interval, and an elementary limit argument. It should be an auxiliary endpoint, not a prerequisite for the extractor theorem. The dyadic iff theorem and fixed-two-bit obstruction are finite counting/parity results and can be formalized independently.

## Completion criteria

Pin a Lean toolchain and Mathlib revision at project creation. Maintain the accompanying ledger as definitions and theorem statements are checked against the manuscript. Record dependencies and any proof changes explicitly. Do not weaken endpoints or omit the auxiliary results while reporting full-paper completion.

A completion report should include a clean full build, no `sorry` or `admit`, no project-owned axioms, and an audit of the axioms used by the public endpoints. Check the actual statements and quantifiers independently of whether Lean accepts their proofs. No such report exists yet for this revision.

## Documentation consulted

Official Mathlib documentation, accessed September 20, 2026:

- [Binomial distributions](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Probability/Distributions/Binomial.html): explicit masses and finite-sum interfaces.
- [Lean tactic reference](https://lean-lang.org/doc/reference/latest/Tactic-Proofs/Tactic-Reference/): proof-producing tactic and computation interfaces.

These links describe existing infrastructure, not a formalization of the present paper. Recheck them against the eventual pinned toolchain.