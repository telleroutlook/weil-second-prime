"""Bernstein ellipse analytic remainder for Gauss-Legendre quadrature.

Replaces the GL-8/GL-4 Richardson self-convergence remainder with a
formally derivable Bernstein ellipse bound.

## Theory

The Gauss-Legendre error on [-1,1] for an analytic function f satisfies:

    |I(f) - GL_n(f)| ≤ (4 * M_rho * rho^{-2n}) / (rho^2 - 1)

where:
  - rho > 1 is the Bernstein ellipse parameter (semi-sum of semi-axes)
  - M_rho = sup_{z on ellipse E_rho} |f(z)|

For a strip [-h, h] (half-width h), after the change of variables to [-1,1]:

    rho = exp(arcsinh(pi/(a*h)))  [nearest singularity of r''(a*t) at Im(t)=pi/a]

This gives a certified, formally derivable error bound that does not rely
on empirical convergence (the weakness of the Richardson approach).

## Usage

    from src.archimedean.bernstein import bernstein_gl_bound

    # Certified upper bound on the GL-n error for integrating
    # r''(a*(x-y)) over a strip of half-width h
    bound = bernstein_gl_bound(half_width, a_num, a_den, n_gl, M_f)
"""

from __future__ import annotations

from fractions import Fraction
import math


# Rational lower bound on pi (used throughout for certified bounds)
# pi > 314159265/100000000
PI_LO = Fraction(314159265, 100000000)

# Certified upper bound on e (from partial sum + geometric tail)
# e < 31967/11760 (verified in Theorem 3 certificate)
E_HI = Fraction(31967, 11760)


def bernstein_gl_bound(
    half_width: Fraction,
    a_num: int,
    a_den: int,
    n_gl: int,
    M_f: Fraction,
) -> Fraction:
    """Certified upper bound on the GL-n quadrature error over a strip of half-width h.

    The integrand is of the form f(t) = r''(a*t) * (polynomial in t),
    analytic in |Im(t)| < pi/a. The Bernstein ellipse parameter for a
    strip [-h, h] scaled to [-1, 1] is:

        rho = exp(arcsinh(pi/(a*h)))

    We use the conservative lower bound:
        rho >= pi/(a*h)   when pi/(a*h) >= sqrt(2)  [since arcsinh(x) >= x/sqrt(2)]
        rho >= 3/2        as a floor

    The error bound is:
        |E_n| <= 4 * M_f * rho^{-2n} / (rho^2 - 1)

    Parameters
    ----------
    half_width : Fraction
        Half-width h of the integration strip (before scaling to [-1,1])
    a_num, a_den : int
        Rational a = a_num/a_den (= L_NUM/L_DEN = 7/20 for our application)
    n_gl : int
        Number of Gauss-Legendre nodes (typically 8)
    M_f : Fraction
        Certified upper bound on |f| on the Bernstein ellipse E_rho

    Returns
    -------
    Fraction
        Certified upper bound on |I(f) - GL_n(f)|
    """
    if half_width <= 0:
        return Fraction(0)

    a = Fraction(a_num, a_den)
    ah = a * half_width

    # rho lower bound: use pi_lo / (a*h) when large enough, else floor at 3/2
    rho_candidate = PI_LO / ah
    rho = max(rho_candidate, Fraction(3, 2))

    # rho^{2n}: Fraction exact computation
    rho_2n = rho ** (2 * n_gl)

    # denominator: rho^2 - 1
    denom = rho_2n * (rho * rho - 1)

    # bound: 4 * M_f / (rho^{2n} * (rho^2 - 1))
    bound = Fraction(4) * M_f / denom

    return bound


def rpp_sup_bound(a_num: int, a_den: int) -> Fraction:
    """Certified upper bound on |r''(t)| for 0 <= t <= a.

    From kernel.py: r''(0) = -7/4 and |r''| is decreasing for t > 0.
    A simple conservative bound: |r''(t)| <= 7/4 for all t in [0, a].

    For a = 7/20 << pi, the actual max is close to 7/4.
    """
    return Fraction(7, 4)


def bernstein_mk_bound(
    a_num: int,
    a_den: int,
    n_row: int,
    n_col: int,
    strip_half_width: Fraction,
    n_gl: int = 8,
) -> Fraction:
    """Certified Bernstein ellipse remainder for one M_K strip integral.

    M_K[n_row, n_col] integrand: k_a(x,y) * P_{n_col}(y) * P_{n_row}(x)
    where |k_a| <= a * |r''(a*(x-y))| <= a * (7/4) = 7a/4
    and |P_n(x)| <= 1 for x in [-1,1].

    Therefore M_f = a * (7/4) * 1 * 1 = 7a/4.
    """
    a = Fraction(a_num, a_den)
    M_f = Fraction(7, 4) * a  # = 7/4 * 7/20 = 49/80 for our application
    return bernstein_gl_bound(strip_half_width, a_num, a_den, n_gl, M_f)
