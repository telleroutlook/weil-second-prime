"""
Directed rational-endpoint interval arithmetic.

All operations use outward rounding:
  - Python Fraction is exact, so arithmetic on Fraction endpoints is exact.
  - When converting FROM Arb balls (float/mpmath), use arb_to_interval() in
    src/constants/constants.py which adds an explicit outward margin.

Reference: proof/schur-certificate.tex, Section 4.
"""

from fractions import Fraction
from typing import Tuple

Interval = Tuple[Fraction, Fraction]


def add(a: Interval, b: Interval) -> Interval:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: Interval, b: Interval) -> Interval:
    return (a[0] - b[1], a[1] - b[0])


def neg(a: Interval) -> Interval:
    return (-a[1], -a[0])


def mul(a: Interval, b: Interval) -> Interval:
    products = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    return (min(products), max(products))


def scalar_mul(c: Fraction, a: Interval) -> Interval:
    if c >= 0:
        return (c * a[0], c * a[1])
    else:
        return (c * a[1], c * a[0])


def div_outward(a: Interval, b: Interval) -> Interval:
    """
    Outward-rounded division a / b.
    Requires b strictly positive (b[0] > 0).
    """
    if b[0] <= 0:
        raise ValueError(f"Divisor interval must be strictly positive, got [{b[0]}, {b[1]}]")
    products = [a[0]/b[0], a[0]/b[1], a[1]/b[0], a[1]/b[1]]
    return (min(products), max(products))


def point(x: Fraction) -> Interval:
    """Create a point interval [x, x]."""
    return (x, x)


def hull(a: Interval, b: Interval) -> Interval:
    """Interval hull (union enclosure)."""
    return (min(a[0], b[0]), max(a[1], b[1]))


def intersect(a: Interval, b: Interval) -> Interval:
    """
    Interval intersection. Raises ValueError if empty.
    Core operation for the Path A / Path B crosscheck.
    """
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    if lo > hi:
        raise ValueError(
            f"Empty intersection: [{float(a[0]):.6e}, {float(a[1]):.6e}] "
            f"∩ [{float(b[0]):.6e}, {float(b[1]):.6e}]"
        )
    return (lo, hi)


def contains(a: Interval, x: Fraction) -> bool:
    """Return True if x is in the closed interval a."""
    return a[0] <= x <= a[1]


def is_strictly_positive(a: Interval) -> bool:
    """Return True iff lower endpoint > 0."""
    return a[0] > 0


def is_positive(a: Interval) -> bool:
    """Return True iff lower endpoint >= 0."""
    return a[0] >= 0


def width(a: Interval) -> Fraction:
    return a[1] - a[0]


def lower(a: Interval) -> Fraction:
    return a[0]


def upper(a: Interval) -> Fraction:
    return a[1]


def sq_lower(a: Interval) -> Fraction:
    """
    Lower bound on x^2 for x in the interval a = [lo, hi].

    If the interval contains zero (lo <= 0 <= hi): the minimum of x^2 is 0.
    Otherwise: the minimum is at the endpoint nearer to zero.

    This is the correct lower bound for the Bessel residual partial norm
    computation. Using lo^2 alone is wrong when lo < 0 < hi.
    """
    lo, hi = a
    if lo <= 0 <= hi:
        return Fraction(0)
    return min(lo * lo, hi * hi)


def sum_intervals(intervals: list) -> Interval:
    """Sum a list of intervals with outward rounding."""
    lo = Fraction(0)
    hi = Fraction(0)
    for iv in intervals:
        lo += iv[0]
        hi += iv[1]
    return (lo, hi)


def mat_vec_mul(M: list, v: list) -> list:
    """
    Multiply a matrix M (list of rows, each row a list of Intervals)
    by a vector v (list of Intervals). Returns a list of Intervals.
    """
    n = len(M)
    result = []
    for i in range(n):
        row_sum = point(Fraction(0))
        for j in range(len(v)):
            row_sum = add(row_sum, mul(M[i][j], v[j]))
        result.append(row_sum)
    return result
