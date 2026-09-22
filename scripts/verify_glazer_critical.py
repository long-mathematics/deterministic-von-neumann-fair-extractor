#!/usr/bin/env python3
"""Exact certificates for the logarithmically improved Glazer extractor.

All mathematical comparisons use integers or fractions. No numerical sampling
or floating-point comparison is used to establish fairness. Default verification
covers the six blocks below the analytic threshold N=128. Larger blocks may be
requested with --extra. The resulting certificate is NOT a Lean formalization.

Examples:
    python verify_glazer_critical.py
    python verify_glazer_critical.py --extra 128 256 512 1024
    python verify_glazer_critical.py --export certificates/finite_blocks.json
    python verify_glazer_critical.py --check-certificate certificates/finite_blocks.json

For N<=128, identities are expanded coefficient by coefficient. For larger N,
a provably injective radix evaluation checks the entire polynomial identity;
it is an exact integer check, not a modular or probabilistic fingerprint.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from fractions import Fraction
from math import comb
from typing import Dict, List, Sequence, Tuple

SMALL_BLOCKS = (2, 4, 8, 16, 32, 64)


def log_interval(x_lo: Fraction, x_hi: Fraction, terms: int) -> Tuple[Fraction, Fraction]:
    """Rigorous interval for log(x), using 2*atanh((x-1)/(x+1))."""
    if not (0 < x_lo <= x_hi):
        raise ValueError("positive ordered interval required")
    z_lo = (x_lo - 1) / (x_lo + 1)
    z_hi = (x_hi - 1) / (x_hi + 1)
    if not (0 <= z_lo <= z_hi < 1):
        raise ValueError("this implementation expects x >= 1")

    def partial(z: Fraction) -> Fraction:
        total = Fraction(0)
        power = z
        z2 = z * z
        for j in range(terms + 1):
            total += 2 * power / (2 * j + 1)
            power *= z2
        return total

    lower = partial(z_lo)
    # The first omitted denominator is 2*terms+3.  Bound every later
    # denominator below by that value and sum the resulting geometric series.
    remainder = (
        2
        * z_hi ** (2 * terms + 3)
        / ((2 * terms + 3) * (1 - z_hi * z_hi))
    )
    upper = partial(z_hi) + remainder
    return lower, upper


def critical_intervals() -> Tuple[Tuple[Fraction, Fraction], Tuple[Fraction, Fraction]]:
    """Return rigorous intervals for beta=C_*-1 and alpha=2-C_*."""
    scale = 10**18
    sqrt5_lo = Fraction(2236067977499789696, scale)
    sqrt5_hi = Fraction(2236067977499789697, scale)
    assert sqrt5_lo * sqrt5_lo < 5 < sqrt5_hi * sqrt5_hi

    phi_lo = (1 + sqrt5_lo) / 2
    phi_hi = (1 + sqrt5_hi) / 2
    log_phi_lo, log_phi_hi = log_interval(phi_lo, phi_hi, 120)
    log5_lo, log5_hi = log_interval(Fraction(5), Fraction(5), 220)

    beta_lo = 2 * log_phi_lo / log5_hi
    beta_hi = 2 * log_phi_hi / log5_lo
    # A compact independently checkable enclosure printed in the manuscript.
    assert beta_lo > Fraction(5979874356654401, 10**16)
    assert beta_hi < Fraction(5979874356654403, 10**16)

    alpha_lo = 1 - beta_hi
    alpha_hi = 1 - beta_lo
    # Use the certified compact enclosure for subsequent floor computations.
    beta_lo = Fraction(5979874356654401, 10**16)
    beta_hi = Fraction(5979874356654403, 10**16)
    return (beta_lo, beta_hi), (1 - beta_hi, 1 - beta_lo)


BETA_INTERVAL, ALPHA_INTERVAL = critical_intervals()


def fibonacci_pair(n: int) -> Tuple[int,int]:
    """Return (F_n,F_{n+1}) by exact integer fast doubling."""
    if n < 0:
        raise ValueError("nonnegative index required")
    if n == 0:
        return 0,1
    a,b=fibonacci_pair(n//2)
    c=a*(2*b-a)
    d=a*a+b*b
    return (d,c+d) if n%2 else (c,d)


def floor_beta(m: int) -> int:
    """Certified floor(beta*m) from comparisons in Q(sqrt(5)).

    5^r < phi^(2m) < 5^(r+1) is equivalent to r < beta*m < r+1.
    No transcendental approximation is needed for these comparisons.
    """
    if m < 1:
        raise ValueError("positive multiplier required")
    f,g=fibonacci_pair(2*m)
    A,B=2*g-f,f  # 2 phi^(2m) = A + B sqrt(5)

    def exceeds(r: int) -> bool:
        gap=2*5**r-A
        if gap <= 0:
            return True
        difference=5*B*B-gap*gap
        if difference == 0:
            raise AssertionError("impossible equality with sqrt(5)")
        return difference > 0

    lo,hi=0,m  # 1 < phi^(2m) < 5^m
    if not exceeds(lo) or exceeds(hi):
        raise AssertionError("invalid floor bracket")
    while hi-lo > 1:
        mid=(lo+hi)//2
        if exceeds(mid):
            lo=mid
        else:
            hi=mid
    return lo


def floor_alpha(m: int) -> int:
    value=m-1-floor_beta(m)  # alpha*m is irrational
    lo,hi=ALPHA_INTERVAL
    fl=(lo*m).numerator//(lo*m).denominator
    fh=(hi*m).numerator//(hi*m).denominator
    if fl == fh and value != fl:
        raise AssertionError(("algebraic and logarithmic floors disagree",m))
    return value


def ceil_cstar(m: int) -> int:
    # C_* m = 2m-alpha*m, and alpha*m is irrational.
    return 2 * m - floor_alpha(m)


def logarithmic_shift(N: int) -> int:
    if N < 2 or N & (N - 1):
        raise ValueError("N must be a power of two, at least 2")
    exponent, power = 0, 1
    while power * 5 <= N:
        exponent += 1
        power *= 5
    return max(0, exponent - 3)


def profile(N: int) -> Tuple[List[int], List[int]]:
    shift = logarithmic_shift(N)
    n = N - 1
    a: List[int] = []
    R: List[int] = []
    for i in range(N):
        ai = max(0, floor_alpha(N + i) - 2 * i + shift)
        ri = n - i - ai
        if not (0 <= ai <= n - i and 0 <= ri <= n - i):
            raise AssertionError(("invalid profile", N, i, ai, ri))
        a.append(ai)
        R.append(ri)
    for i in range(N - 1):
        if a[i] - a[i + 1] not in (0, 1, 2):
            raise AssertionError(("invalid drop", N, i, a[i], a[i + 1]))
    return a, R


def folded_weights(N: int, R: Sequence[int]) -> List[Dict[int, int]]:
    """Signed numerators of W_i, with common denominator 2^(N-1)."""
    n = N - 1
    positive = [comb(n, t) for t in range(n + 1)]
    sign = 1
    output: List[Dict[int, int]] = []

    for i in range(N):
        if not positive:
            output.append({})
            continue
        degree = max(j for j, value in enumerate(positive) if value)
        Ri = R[i]
        if degree <= Ri:
            output.append(
                {t: sign * value for t, value in enumerate(positive) if value}
            )
            positive = []
            continue

        tail = sum(positive[Ri:])
        digit = {t: sign * positive[t] for t in range(Ri) if positive[t]}
        if tail:
            digit[Ri] = sign * tail
        output.append(digit)

        next_positive = [0] * degree
        suffix = 0
        for s in range(degree - 1, Ri - 1, -1):
            suffix += positive[s + 1]
            next_positive[s] = 2 * suffix
        positive = next_positive
        sign = -sign

    if positive:
        raise AssertionError(("fold did not terminate", N))
    return output


def multiply(a: Sequence[int], b: Sequence[int]) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        if av:
            for j, bv in enumerate(b):
                if bv:
                    out[i + j] += av * bv
    return out


def parity_character_polynomial(R: int, t: int) -> List[int]:
    """Coefficients of (1-x)^t (1+x)^(R-t), by its differential recurrence."""
    if not 0 <= t <= R:
        raise ValueError("require 0 <= t <= R")
    result = [1]
    previous = 0
    for r in range(R):
        numerator = (R - 2*t) * result[-1] - (R-r+1) * previous
        value, remainder = divmod(numerator, r+1)
        if remainder:
            raise AssertionError(("nonintegral character recurrence", R,t,r))
        previous = result[-1]
        result.append(value)
    return result


def real_center(
    N: int, a: Sequence[int], R: Sequence[int]
) -> Tuple[List[List[Fraction]], List[Dict[int, int]]]:
    denominator = 2 ** (N - 1)
    weights = folded_weights(N, R)
    center: List[List[Fraction]] = []
    for i, (Ri, digit) in enumerate(zip(R, weights)):
        e = [0] * (Ri + 1)
        for t, numerator in digit.items():
            character = parity_character_polynomial(Ri, t)
            for r, coefficient in enumerate(character):
                e[r] += numerator * coefficient
        # Substitute x=-z.
        center.append(
            [Fraction((-1) ** (i + r) * coefficient, denominator) for r, coefficient in enumerate(e)]
        )
    return center, weights


def nearest_parity_integer(x: Fraction, parity: int, boundary: bool = False) -> int:
    if boundary:
        candidates = (-1, 1)
    else:
        floor_x = x.numerator // x.denominator
        candidates = tuple(j for j in range(floor_x - 4, floor_x + 5) if j % 2 == parity)
    return min(candidates, key=lambda j: (abs(Fraction(j) - x), abs(j), j))


def integral_rounding(
    N: int,
    a: Sequence[int],
    R: Sequence[int],
    u: Sequence[Sequence[Fraction]],
) -> Tuple[List[List[int]], Fraction]:
    """Construct the exact parity-correct packet coefficients."""
    delta: List[List[Fraction]] = [
        [Fraction(0)] * (R[i] + 1) for i in range(N - 1)
    ]
    f: List[List[int]] = [[0] * (R[i] + 1) for i in range(N)]
    max_error = Fraction(0)

    def predecessor(i: int, r: int) -> Fraction:
        if i == 0:
            return Fraction(0)
        s = a[i - 1] - a[i]
        value = Fraction(0)
        for j in range(s + 1):
            parent_r = r + 1 - j
            if 1 <= parent_r <= R[i - 1]:
                value += (-1) ** j * comb(s, j) * delta[i - 1][parent_r]
        return value

    for i in range(N):
        # Lower packet boundary.  It has already been forced by delta[i-1,1].
        lower = u[i][0] - predecessor(i, 0)
        if lower.denominator != 1 or abs(lower) != 1:
            raise AssertionError(("lower boundary", N, i, lower))
        f[i][0] = int(lower)
