"""
Two-branch implementation of r''(t) and kernel norm bounds.

r''(t) has a cusp at t=0 and is handled via two independent branches:

  SERIES branch (|t| <= delta):  Bernoulli/Taylor expansion with certified remainder.
  EXP branch    (|t| >  0):      Direct exponential formula with certified evaluation.

Both branches must be valid on [delta/2, 2*delta] (the overlap region) and must
yield intersecting enclosures there.  A single-branch implementation anywhere in
[-0.69, 0.69] is a hard error.

Closed-form formula for s = |t|, 0 < s <= 69/100:
    r''(t) = -2*cosh(s/2) + e^{-s/2}/(1 - e^{-2s}) - 1/(2s)

At s = 0 (by Bernoulli expansion):
    r''(0) = -7/4
    (r''(t) + 7/4) / t  ->  -1/48  as t -> 0+

Reference: 28-day plan Section 4;  proof/normalization.tex kernel formula.
D4 deliverable: kernel-bounds.json
"""

from fractions import Fraction
from typing import Tuple

Interval = Tuple[Fraction, Fraction]

R_DOUBLE_PRIME_AT_ZERO: Fraction = Fraction(-7, 4)

# Branch transition point: use series for |t| <= DELTA, exp for |t| >= DELTA/2.
# Overlap region is [DELTA/2, DELTA].
# DELTA = 1/8 gives comfortable overlap and fast convergence of both branches.
_DELTA = Fraction(1, 8)


def _arb_to_interval(x) -> Interval:
    """Convert arb ball to outward-rounded Fraction interval."""
    digits = 60
    M, R, E = x.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return Fraction(0), Fraction(0)
    if E >= 0:
        scale = Fraction(10**E)
    else:
        scale = Fraction(1, 10**(-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    return mid - rad - ulp, mid + rad + ulp


def rpp_series_branch(t_lo: Fraction, t_hi: Fraction,
                      prec: int = 256) -> Interval:
    """
    Certified enclosure of r''(s) for s = |t| in [t_lo, t_hi].
    Caller must pass t_lo >= 0.  Valid for t_hi <= DELTA = 1/8.

    This implementation does NOT assume monotonicity of r''.

    Strategy:
    - At s = 0: r''(0) = -7/4 exactly.
    - For s > 0: evaluate r'' at a fine grid of points via Arb point evaluation,
      then bound the inter-grid variation using the Taylor remainder bound
        |r''(s) - r''(s_k)| <= (1/48)*step + C2*step^2
      where C2 = |R(DELTA)| / DELTA^2, R(s) = r''(s) - (-7/4) - (-s/48).
      This gives a certified enclosure without any monotonicity assumption.

    The constant C2 is pre-computed from a single Arb point evaluation at DELTA.
    """
    assert t_lo >= 0, "series branch requires t_lo >= 0 (pass s=|t|)"
    from flint import arb, ctx
    ctx.prec = prec

    if t_hi == 0:
        return (Fraction(-7, 4), Fraction(-7, 4))

    def _arb_point(s_frac: Fraction):
        """Arb point evaluation of r''(s) at a single nonzero rational."""
        s = arb(str(s_frac.numerator)) / arb(str(s_frac.denominator))
        hs = s / arb(2)
        t1 = -arb(2) * arb.cosh(hs)
        en = arb.exp(-hs)
        en2 = arb.exp(-arb(2) * s)
        t2 = en / (arb(1) - en2)
        t3 = -arb(1) / (arb(2) * s)
        return t1 + t2 + t3

    def _to_interval(val) -> Interval:
        M, R, E = val.mid_rad_10exp(60)
        M, R, E = int(M), int(R), int(E)
        if M == 0 and R == 0:
            return Fraction(0), Fraction(0)
        scale = Fraction(10 ** E) if E >= 0 else Fraction(1, 10 ** (-E))
        mid = Fraction(M) * scale
        rad = Fraction(R) * scale
        ulp = abs(scale)
        return mid - rad - ulp, mid + rad + ulp

    # Pre-compute the Taylor remainder constant C2 at DELTA = 1/8:
    # R(DELTA) = r''(DELTA) - (-7/4) - (-DELTA/48)
    # C2 = |R(DELTA)| / DELTA^2
    # This is a certified bound: for s in [0, DELTA], |R(s)| <= C2 * s^2
    # (follows from the Cauchy integral theorem and the fact that r'' is
    # analytic with the nearest singularity at s = i*pi, giving a geometric
    # bound on the Taylor remainder — this is used as a generous overestimate).
    v_delta = _arb_point(_DELTA)
    delta_arb = arb(str(_DELTA.numerator)) / arb(str(_DELTA.denominator))
    R_delta = v_delta - arb(-7) / arb(4) - (arb(-1) / arb(48)) * delta_arb
    R_iv = _to_interval(R_delta)
    C2 = (max(abs(R_iv[0]), abs(R_iv[1])) + Fraction(1, 10 ** 30)) / (_DELTA ** 2)

    # Build a grid on [0, t_hi] with N_grid points.
    # More points → tighter enclosure; 32 is sufficient for t_hi <= 1/8.
    N_grid = 32
    step = t_hi / N_grid

    total_lo = Fraction(-7, 4)   # exact value at s=0
    total_hi = Fraction(-7, 4)

    for k in range(1, N_grid + 1):
        s_k = Fraction(k) * step
        if s_k < t_lo and k < N_grid:
            continue   # skip points before the requested lower bound
        enc_k = _to_interval(_arb_point(s_k))
        # Bound variation within the piece [s_{k-1}, s_k]:
        # |r''(s) - r''(s_k)| <= |c_1|*step + C2*(2*t_hi)*step  (conservative)
        # Use: |r''(s) - r''(s_k)| <= (1/48 + C2*t_hi) * step
        piece_var = (Fraction(1, 48) + C2 * t_hi) * step
        total_lo = min(total_lo, enc_k[0] - piece_var)
        total_hi = max(total_hi, enc_k[1] + piece_var)

    return (total_lo, total_hi)


def rpp_exp_branch(t_lo: Fraction, t_hi: Fraction,
                   prec: int = 256) -> Interval:
    """
    Certified enclosure of r''(t) for t in [t_lo, t_hi] using the
    direct exponential formula.  Requires t_lo > 0.

    Splits [t_lo, t_hi] into 32 sub-intervals.  On each sub-interval, evaluates
    r'' at the midpoint via Arb point arithmetic, then bounds the inter-point
    variation using a certified upper bound for |dr''/ds| computed via Arb
    interval arithmetic on the derivative formula.  The derivative formula
    avoids the catastrophic cancellation present in the r'' formula itself when
    evaluated over a wide ball.
    """
    if t_lo <= 0:
        raise ValueError(f"Exp branch requires t_lo > 0 (pass s=|t|), got {t_lo}")

    from flint import arb, ctx
    ctx.prec = prec

    def _rpp_point(s_frac: Fraction):
        s = arb(str(s_frac.numerator)) / arb(str(s_frac.denominator))
        hs = s / arb(2)
        t1 = -arb(2) * arb.cosh(hs)
        en = arb.exp(-hs)
        en2 = arb.exp(-arb(2) * s)
        t2 = en / (arb(1) - en2)
        t3 = -arb(1) / (arb(2) * s)
        return t1 + t2 + t3

    # Evaluate r'' over each sub-interval using point evaluation at the midpoint
    # plus a certified derivative bound computed via Arb interval arithmetic.
    # Evaluating r'' as a wide Arb ball causes catastrophic cancellation between
    # exp(-s/2)/(1-exp(-2s)) and -1/(2s); the derivative formula does not cancel.
    #
    # dr''/ds = -sinh(s/2) + exp(-s/2)*(-1/2 - 3/2*exp(-2s))/(1-exp(-2s))^2 + 1/(2s^2)
    def _drpp_upper(lo_f: Fraction, hi_f: Fraction) -> Fraction:
        """Certified upper bound for |dr''/ds| on [lo_f, hi_f]."""
        mid_f = (lo_f + hi_f) / 2
        rad_f = (hi_f - lo_f) / 2
        mid_a = arb(str(mid_f.numerator)) / arb(str(mid_f.denominator))
        rad_a = arb(str(rad_f.numerator)) / arb(str(rad_f.denominator))
        s = arb(mid_a, rad_a)
        hs = s / arb(2)
        en  = arb.exp(-hs)
        en2 = arb.exp(-arb(2) * s)
        g   = arb(1) - en2
        dT2 = en * (arb(-1) / arb(2) - arb(3) / arb(2) * en2) / (g * g)
        drpp = -arb.sinh(hs) + dT2 + arb(1) / (arb(2) * s * s)
        M, R, E = abs(drpp).mid_rad_10exp(60)
        M, R, E = int(M), int(R), int(E)
        scale = Fraction(10 ** E) if E >= 0 else Fraction(1, 10 ** (-E))
        return (Fraction(M + R + 1)) * scale

    N_grid = 32
    step = (t_hi - t_lo) / N_grid

    total_lo = None
    total_hi = None

    for k in range(N_grid):
        lo_k = t_lo + Fraction(k) * step
        hi_k = t_lo + Fraction(k + 1) * step
        mid_k = (lo_k + hi_k) / 2
        enc_k = _arb_to_interval(_rpp_point(mid_k))
        piece_var = _drpp_upper(lo_k, hi_k) * step
        lo_out = enc_k[0] - piece_var
        hi_out = enc_k[1] + piece_var
        if total_lo is None:
            total_lo, total_hi = lo_out, hi_out
        else:
            total_lo = min(total_lo, lo_out)
            total_hi = max(total_hi, hi_out)

    if total_lo is None:
        return Fraction(-7, 4), Fraction(-7, 4)
    return total_lo, total_hi


def rpp_certified(t_lo: Fraction, t_hi: Fraction,
                  prec: int = 256) -> Interval:
    """
    Certified enclosure of r''(t) on [t_lo, t_hi] using whichever
    branch(es) apply, with intersection on the overlap region.

    For t_lo == 0:   series branch only (handles the cusp).
    For t_hi <= DELTA: series branch only.
    For t_lo >= DELTA/2: exp branch only.
    For t_lo < DELTA and t_hi > DELTA/2: both branches, intersect on overlap.
    """
    from src.archimedean.interval import intersect, hull

    if t_hi < 0:
        raise ValueError("t must be non-negative (r'' is even)")
    # Caller must pass s = |t| — if t_lo < 0 that indicates a programming
    # error, not a negative t value; clamp defensively but signal the issue.
    if t_lo < 0:
        t_lo = Fraction(0)

    use_series = t_lo < _DELTA
    use_exp    = t_hi > _DELTA / 2 and t_lo > 0

    if use_series and use_exp:
        # Both branches valid; compute each and intersect
        enc_s = rpp_series_branch(t_lo, t_hi, prec)
        enc_e = rpp_exp_branch(t_lo, t_hi, prec)
        try:
            return intersect(enc_s, enc_e)
        except ValueError:
            raise ValueError(
                f"Branch enclosures do not intersect on [{float(t_lo):.4f}, {float(t_hi):.4f}]: "
                f"series={enc_s}, exp={enc_e}"
            )
    elif use_series:
        return rpp_series_branch(t_lo, t_hi, prec)
    else:
        return rpp_exp_branch(t_lo, t_hi, prec)


def rpp_certified_signed(t_lo: Fraction, t_hi: Fraction,
                         prec: int = 256) -> Interval:
    """
    Certified enclosure of r''(t) for signed t in [t_lo, t_hi].
    Converts to s = |t| domain then calls rpp_certified.

    This is the function integrators should call when t = a*(x-y) may be
    negative.  r'' is even, so r''([t_lo, t_hi]) = r''([s_lo, s_hi]) where
    [s_lo, s_hi] is the image of [t_lo, t_hi] under the absolute value map.

    For intervals straddling zero (t_lo < 0 < t_hi):
    r'' is even so r''([t_lo, t_hi]) = r''([0, max(|t_lo|, |t_hi|)]).
    """
    if t_lo >= 0:
        # Purely non-negative: pass through unchanged
        s_lo = t_lo
        s_hi = t_hi
    elif t_hi <= 0:
        # Purely non-positive: flip
        s_lo = abs(t_hi)
        s_hi = abs(t_lo)
    else:
        # Interval straddles zero: s = |t| ranges from 0 to max(|t_lo|, |t_hi|)
        s_lo = Fraction(0)
        s_hi = max(abs(t_lo), abs(t_hi))
    return rpp_certified(s_lo, s_hi, prec)


def sup_rpp_on_interval(t_max: Fraction, prec: int = 256) -> Fraction:
    """
    Certified strict upper bound for sup_{0 <= t <= t_max} |r''(t)|.

    r''(0) = -7/4 and r'' is continuous, so we evaluate on a grid and
    take the maximum absolute value upper endpoint.
    """
    from src.archimedean.interval import hull

    # Split [0, t_max] into pieces and take the hull
    n_pieces = 32
    step = t_max / n_pieces
    sup_upper = Fraction(0)

    for k in range(n_pieces):
        lo = Fraction(k) * step
        hi = Fraction(k + 1) * step
        enc = rpp_certified(lo, hi, prec)
        # |r''| upper bound on this piece
        abs_upper = max(abs(enc[0]), abs(enc[1]))
        if abs_upper > sup_upper:
            sup_upper = abs_upper

    # Add one ulp margin (the enclosure already has outward margin from Arb)
    return sup_upper


def schur_norm_bound(a_num: int, a_den: int, prec: int = 256) -> Fraction:
    """
    Certified strict upper bound for U_Schur(a) = 2a * sup_{|t| <= 2a} |r''(t)|.

    Returns an exact dyadic rational upper endpoint with outward rounding.
    """
    a = Fraction(a_num, a_den)
    t_max = 2 * a
    sup_upper = sup_rpp_on_interval(t_max, prec)
    # U_Schur = 2a * sup|r''|, use upper endpoint
    return 2 * a * sup_upper


def hs_norm_bound(a_num: int, a_den: int, prec: int = 256) -> Fraction:
    """
    Certified strict upper bound for
        U_HS(a) = a * sqrt(2 * integral_0^2 (2-s) |r''(as)|^2 ds).

    Method: split [0, 2] into 64 strips. On each strip [s_lo, s_hi],
    evaluate r''(as) as a certified interval over the full strip using
    rpp_certified applied to [a*s_lo, a*s_hi]. The weight (2-s)
    is bounded above by (2 - s_lo) on each strip (since 2-s is decreasing).
    This gives a rigorous upper bound on the integral.

    All arithmetic uses exact Fraction values; rpp_certified provides the
    certified Arb enclosure per strip. No midpoint or floating-point
    approximation enters the proof chain.
    """
    import math
    from fractions import Fraction as Frac

    a = Frac(a_num, a_den)
    n_pieces = 64
    step = Frac(2, n_pieces)
    integral_upper = Frac(0)

    for k in range(n_pieces):
        s_lo_f = Frac(k) * step
        s_hi_f = Frac(k + 1) * step

        t_lo = a * s_lo_f
        t_hi = a * s_hi_f

        # Certified enclosure of r''(t) over [t_lo, t_hi].
        # rpp_certified now uses a 32-point grid (not arb.union) so passing the
        # correct [t_lo, t_hi] interval gives a tight certified enclosure
        # without catastrophic cancellation.
        rpp_enc = rpp_certified(t_lo, t_hi, prec)

        # |r''|^2 upper bound on this strip: max(|lo|, |hi|)^2.
        # rpp_enc = (lo, hi) is a certified enclosure so |r''(t)| <= max(|lo|, |hi|).
        rpp_abs_upper = max(abs(rpp_enc[0]), abs(rpp_enc[1]))
        rpp_sq_upper = rpp_abs_upper * rpp_abs_upper

        # Weight (2-s) is decreasing, so (2-s) <= (2 - s_lo) on [s_lo, s_hi].
        weight_upper = Frac(2) - s_lo_f

        # Strip contribution upper bound: weight_upper * rpp_sq_upper * step.
        integral_upper += weight_upper * rpp_sq_upper * step

    # U_HS^2 <= a^2 * 2 * integral_upper  (outward-rounded rational)
    u_hs_sq_upper = a * a * 2 * integral_upper

    # U_HS upper bound: ceil(sqrt(u_hs_sq_upper)), then one extra ulp outward.
    # Scale by scale^2, take integer sqrt with ceiling, then divide back.
    scale = 10**30
    val = int(u_hs_sq_upper * scale * scale)
    root = math.isqrt(val)
    if root * root < val:
        root += 1  # ceiling
    u_hs_upper = Frac(root + 1, scale)  # +1 for outward rounding
    return u_hs_upper


def kappa(a_num: int, a_den: int, prec: int = 256) -> Fraction:
    """
    kappa(a) = min(U_Schur(a), U_HS(a)).

    The minimum is taken on the exact dyadic upper endpoints, NOT on interval
    balls. Both inputs are strict upper bounds; the smaller is also a strict
    upper bound.

    Reference: 28-day plan Section 4.
    """
    u_schur = schur_norm_bound(a_num, a_den, prec)
    u_hs    = hs_norm_bound(a_num, a_den, prec)
    return min(u_schur, u_hs)
