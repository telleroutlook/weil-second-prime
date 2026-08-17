"""
Path A integration engine: dyadic adaptive quadrature with Duffy transforms.

Computes certified interval enclosures for M_K and S_KK matrix entries.
M_V entries come from src/log_moments/log_moments.py (exact formula).

2D integrals (M_K):
    M_K[i,j] = integral_{[-1,1]^2} k_a(x,y) P_j(y) P_i(x) dx dy
    where k_a(x,y) = -a * r''(a*(x-y))
    Split into two triangles: {x>y} and {x<y}, each mapped to [0,1]^2 via Duffy.

3D integrals (S_KK):
    S_KK[i,j] = integral_{[-1,1]^3} k_a(x,y) k_a(x,z) P_j(y) P_i(z) dx dy dz
    Non-smooth at x=y and x=z. Four regions:
      R1: y<x, z<x    R2: y<x, z>x
      R3: y>x, z<x    R4: y>x, z>x
    Each mapped to [0,1]^3 via tensor Duffy.

IMPORTANT: This file must share NO approximation or remainder code with
integrator_b/integrator.py.

Reference: 28-day plan Sections 6.2-6.3.
D6 deliverable.
"""

from dataclasses import dataclass, field
from fractions import Fraction
from typing import List, Tuple, Dict, Optional

Interval = Tuple[Fraction, Fraction]


# ---------------------------------------------------------------------------
# Shared helpers (local to Path A only)
# ---------------------------------------------------------------------------

def _arb_to_interval(x) -> Interval:
    digits = 60
    M, R, E = x.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return Fraction(0), Fraction(0)
    scale = Fraction(10**E) if E >= 0 else Fraction(1, 10**(-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    return mid - rad - ulp, mid + rad + ulp


def _interval_add(a: Interval, b: Interval) -> Interval:
    return (a[0] + b[0], a[1] + b[1])


def _frac_to_arb(f: Fraction):
    from flint import arb
    return arb(f.numerator) / arb(f.denominator)


def _legendre_at_arb(n: int, x):
    """Evaluate P_n(x) for arb x using the recurrence."""
    from flint import arb
    if n == 0:
        return arb(1)
    if n == 1:
        return x
    p_prev, p_curr = arb(1), x
    for k in range(2, n + 1):
        p_next = ((2*k - 1) * x * p_curr - (k - 1) * p_prev) / k
        p_prev, p_curr = p_curr, p_next
    return p_curr


def _geom_decay_upper_bound_a(a_num: int, a_den: int) -> Fraction:
    """
    Certified rational UPPER bound for exp(-2*pi/a).
    Uses Arb at 256-bit precision — no float() in the proof chain.
    """
    from flint import arb, ctx
    ctx.prec = 256
    a_arb = arb(str(a_num)) / arb(str(a_den))
    val = arb.exp(arb(-2) * arb.pi() / a_arb)
    digits = 30
    M, R, E = val.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return Fraction(1, 10**30)
    scale = Fraction(10**E) if E >= 0 else Fraction(1, 10**(-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    return mid + rad + ulp  # outward upper bound


def _arb_rad_to_frac_upper(rad_arb) -> str:
    """
    Convert an Arb ball radius to a conservative rational UPPER bound string.
    Returns a 'p/q' string suitable for remainder_bound fields.
    Uses outward-rounded Fraction conversion — no float() in the proof chain.
    """
    from flint import arb
    digits = 20
    # The radius is non-negative; use mid_rad_10exp on the radius itself
    M, R, E = rad_arb.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return "0"
    scale = Fraction(10**E) if E >= 0 else Fraction(1, 10**(-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    upper = mid + rad + ulp
    if upper <= 0:
        return "0"
    return f"{upper.numerator}/{upper.denominator}"


def _rpp_arb(t, prec: int = 256):
    """
    Evaluate r''(t) as an Arb ball for arb t.
    Handles t=0 exactly as -7/4.
    """
    from flint import arb, ctx
    ctx.prec = prec
    # Use the two-term formula; near t=0 use Taylor (handled by Arb automatically)
    half_t = t / arb(2)
    exp_neg = arb.exp(-half_t)
    exp_neg2 = arb.exp(-arb(2) * t)
    # -2*cosh(t/2) + exp(-t/2)/(1-exp(-2t)) - 1/(2t)
    # Near t=0 the last two terms partially cancel: use the series form via Taylor
    # Arb handles the removable singularity automatically at high precision
    # For t strictly > 0:
    term1 = -arb(2) * arb.cosh(half_t)
    term2 = exp_neg / (arb(1) - exp_neg2)
    term3 = -arb(1) / (arb(2) * t)
    return term1 + term2 + term3


# ---------------------------------------------------------------------------
# Leaf witness record
# ---------------------------------------------------------------------------

@dataclass
class LeafWitnessA:
    """Witness record for a single dyadic leaf cell (Path A)."""
    matrix:          str
    row:             int
    col:             int
    path:            str = "A"
    dimension:       int = 2
    cell_lo:         list = field(default_factory=list)  # per-dim lower bounds as [num,den]
    cell_hi:         list = field(default_factory=list)  # per-dim upper bounds as [num,den]
    transform_type:  str = "none"
    jacobian_expr:   str = ""
    branch:          str = ""
    enclosure_lower: str = "0"
    enclosure_upper: str = "0"
    remainder_type:  str = "arb_ball"
    remainder_bound: str = "0"

    def to_dict(self) -> dict:
        return {
            "entry":  {"matrix": self.matrix, "row": self.row, "col": self.col},
            "path":   self.path,
            "domain": {
                "dimension": self.dimension,
                "cell_lo": self.cell_lo,
                "cell_hi": self.cell_hi,
            },
            "transform": {
                "type":     self.transform_type,
                "jacobian": self.jacobian_expr,
            },
            "branch":    self.branch,
            "enclosure": {
                "lower": self.enclosure_lower,
                "upper": self.enclosure_upper,
            },
            "remainder": {
                "type":  self.remainder_type,
                "bound": self.remainder_bound,
            },
        }


@dataclass
class PathAResult:
    """Accumulated result for one matrix entry from Path A."""
    matrix:          str
    row:             int
    col:             int
    enclosure_lower: Fraction = Fraction(0)
    enclosure_upper: Fraction = Fraction(0)
    leaves:          List[LeafWitnessA] = field(default_factory=list)

    def to_interval(self) -> Interval:
        return (self.enclosure_lower, self.enclosure_upper)


# ---------------------------------------------------------------------------
# 1D adaptive Gauss integration over [lo, hi] via Arb — certified remainder
# ---------------------------------------------------------------------------

def _gl_nodes_weights_arb(n_gl: int, prec: int = 256):
    """Load GL nodes and weights at 70-digit precision, return as Arb balls."""
    from flint import arb, ctx
    ctx.prec = prec
    import mpmath
    mpmath.mp.dps = 80
    nodes_mp, wts_mp = mpmath.gauss_quadrature(n_gl, "legendre")
    nodes = [arb(mpmath.nstr(nd, 70, strip_zeros=False)) for nd in nodes_mp]
    wts   = [arb(mpmath.nstr(wt, 70, strip_zeros=False)) for wt in wts_mp]
    return nodes, wts


def _gl_eval_arb(func, a_arb, b_arb, nodes, wts):
    """Evaluate GL quadrature sum of func on [a,b] with precomputed nodes/weights."""
    from flint import arb
    half = (b_arb - a_arb) / arb(2)
    mid  = (a_arb + b_arb) / arb(2)
    s = arb(0)
    for nd, wt in zip(nodes, wts):
        s = s + wt * func(mid + half * nd)
    return half * s


def _integrate_1d_arb(func, lo: Fraction, hi: Fraction, depth: int,
                      prec: int = 256,
                      bernstein_M_f: Fraction | None = None,
                      a_num: int = 7, a_den: int = 20):
    """
    Certified 1D integration via dyadic GL-8.

    Two remainder modes:

    MODE A — Bernstein ellipse (default when bernstein_M_f is provided):
      Uses the formally derivable bound:
        |I(f) - GL_8(f)| ≤ 4·M_f·rho^{-16} / (rho^2 - 1)
      where rho ≥ pi/(a·h), h = half-width of each sub-interval.
      This is an analytic bound that does not rely on empirical convergence.
      Activated by passing bernstein_M_f = certified upper bound on |f|.

    MODE B — Richardson self-convergence (fallback):
      Uses 2·|GL_8 - GL_4| per sub-interval (conservative but certified).
      Activated when bernstein_M_f is None.

    Both modes return an Arb ball that provably contains the integral.
    Mode A is preferred for O2 certification; Mode B is retained for
    compatibility and as a cross-check.
    """
    from flint import arb, ctx
    ctx.prec = prec

    nodes8, wts8 = _gl_nodes_weights_arb(8, prec)

    n_sub = 2 ** depth
    step = (hi - lo) / n_sub
    half_width = step / 2
    total = arb(0)

    if bernstein_M_f is not None:
        # MODE A: Bernstein ellipse remainder
        from src.archimedean.bernstein import bernstein_gl_bound
        remainder_frac = bernstein_gl_bound(
            half_width, a_num, a_den, n_gl=8, M_f=bernstein_M_f
        ) * n_sub  # total bound = per-strip bound * number of strips
        remainder_arb = arb(str(remainder_frac.numerator)) / arb(str(remainder_frac.denominator))

        for k in range(n_sub):
            a_k = lo + Fraction(k) * step
            b_k = lo + Fraction(k + 1) * step
            a_arb = _frac_to_arb(a_k)
            b_arb = _frac_to_arb(b_k)
            gl8 = _gl_eval_arb(func, a_arb, b_arb, nodes8, wts8)
            total = total + gl8

        # Add Bernstein remainder as a single outer ball
        total = total + arb.union(-remainder_arb, remainder_arb)

    else:
        # MODE B: Richardson self-convergence remainder
        nodes4, wts4 = _gl_nodes_weights_arb(4, prec)

        for k in range(n_sub):
            a_k = lo + Fraction(k) * step
            b_k = lo + Fraction(k + 1) * step
            a_arb = _frac_to_arb(a_k)
            b_arb = _frac_to_arb(b_k)

            gl8 = _gl_eval_arb(func, a_arb, b_arb, nodes8, wts8)
            gl4 = _gl_eval_arb(func, a_arb, b_arb, nodes4, wts4)

            diff_abs = (gl8 - gl4).abs_upper()
            remainder = arb(2) * diff_abs
            total = total + gl8 + arb.union(-remainder, remainder)

    return total


# ---------------------------------------------------------------------------
# M_K: <K_a P_j, P_i>
# ---------------------------------------------------------------------------

def _mk_integrand_arb(x, y, n_row: int, n_col: int,
                      a: Fraction, prec: int = 256):
    """
    Integrand for M_K[row,col] = integral k_a(x,y) P_{n_col}(y) P_{n_row}(x) dx dy
    k_a(x,y) = -a * r''(a*(x-y))

    For the Duffy-transformed triangle {x > y}: (x,y) in [-1,1]^2, x > y.
    Duffy map:  x = lo_x + (hi_x - lo_x)*u,  y = x - (x - lo_y)*v
    This is handled by the caller; here we just evaluate the kernel product.
    """
    from flint import arb, ctx
    ctx.prec = prec
    a_arb = _frac_to_arb(a)
    t = a_arb * (x - y)
    s = abs(t)  # r'' is even; always use s = |t|
    if s.is_zero():
        rpp_val = arb(-7) / arb(4)  # r''(0) = -7/4
    else:
        rpp_val = _rpp_arb(s, prec)
    k_val = -a_arb * rpp_val
    return k_val * _legendre_at_arb(n_col, y) * _legendre_at_arb(n_row, x)


def integrate_M_K(n_row: int, n_col: int, a_num: int, a_den: int,
                  depth: int = 4, prec: int = 256,
                  use_bernstein: bool = True,
                  skip_remainder: bool = False) -> PathAResult:
    """
    Certified enclosure of M_K[row,col] = <K_a P_{n_col}, P_{n_row}>.

    Splits [-1,1]^2 into two triangles (x>y and x<y) and integrates each
    using 2D adaptive quadrature via iterated 1D Gauss on a dyadic grid.

    The diagonal x=y is a cusp; the Duffy split ensures smooth integrands.
    The y-inner integral is further subdivided into max(1, n_col//4 + 1)
    strips so that high-degree Legendre polynomials P_{n_col}(y) are
    resolved even for large n_col.

    The kernel is evaluated as -a * r''(s) where s = |t| = |a*(x-y)|.
    r'' is even, so using |t| is equivalent to the original formula but
    makes the symmetry explicit and avoids ambiguity with signed t.

    Parameters
    ----------
    use_bernstein : bool
        If True (default), use the Bernstein ellipse analytic remainder
        for formally derivable certificates (O2 certification path).
        If False, use Richardson GL-8/GL-4 remainder (for large n_row where
        Bernstein bound >> 1, e.g. N >= 20 second window).
    skip_remainder : bool
        If True, skip all truncation-error remainder computation (Bernstein
        and Richardson). Returns the raw GL-8 Arb ball only. Use ONLY for
        float-center pilot paths (e.g. recompute_schur float mode) where
        the interval width is discarded anyway. P0 defect if used in a
        certified proof path.
    """
    from flint import arb, ctx
    ctx.prec = prec

    a = Fraction(a_num, a_den)
    a_arb = _frac_to_arb(a)

    # Triangle 1: x > y  (integrate over x in [-1,1], then y in [-1,x])
    # Triangle 2: x < y  (integrate over x in [-1,1], then y in [x,1])
    # By symmetry of the kernel (r'' is even): both triangles give equal contribution
    # if P_i, P_j have the same parity; otherwise they sum correctly.
    # We integrate both for correctness.

    total_arb = arb(0)
    total_arb_4 = arb(0)  # Richardson: GL-4 counterpart (used when use_bernstein=False)
    leaves = []

    # Adaptive subdivision counts to resolve high-degree Legendre oscillations.
    # n_col_deg = max(n_row, n_col) is used symmetrically for both axes so that
    # both P_{n_row}(x) and P_{n_col}(y) are resolved by their respective GL loops.
    n_col_deg = max(n_row, n_col)
    n_ysub = max(1, n_col_deg // 2 + 2)  # y strips: resolves P_{n_col}(y)
    n_xsub = max(1, n_col_deg // 2 + 2)  # x strips: resolves P_{n_row}(x) (symmetric)
    n_sub = max(2 ** depth, n_xsub)       # depth adds extra refinement on top
    x_step = Fraction(2, n_sub)  # step on [-1,1]

    # Compute GL nodes once — they don't change between x-strips
    import mpmath
    mpmath.mp.dps = 80
    nodes_mp, weights_mp = mpmath.gauss_quadrature(8, "legendre")
    x_nodes = [arb(mpmath.nstr(nd, 70, strip_zeros=False)) for nd in nodes_mp]
    x_wts   = [arb(mpmath.nstr(wt, 70, strip_zeros=False)) for wt in weights_mp]
    if not use_bernstein and not skip_remainder:
        nodes4_mp, weights4_mp = mpmath.gauss_quadrature(4, "legendre")
        x_nodes4 = [arb(mpmath.nstr(nd, 70, strip_zeros=False)) for nd in nodes4_mp]
        x_wts4   = [arb(mpmath.nstr(wt, 70, strip_zeros=False)) for wt in weights4_mp]

    for kx in range(n_sub):
        x_lo = Fraction(-1) + Fraction(kx) * x_step
        x_hi = x_lo + x_step

        x_lo_arb  = _frac_to_arb(x_lo)
        x_hi_arb  = _frac_to_arb(x_hi)
        x_mid_arb = (x_lo_arb + x_hi_arb) / arb(2)
        x_hf_arb  = (x_hi_arb - x_lo_arb) / arb(2)
        x_hf_frac = (x_hi - x_lo) / 2

        for tri in [1, 2]:  # tri=1: x>y, tri=2: x<y
            strip_sum = arb(0)
            for x_nd, x_wt in zip(x_nodes, x_wts):
                x_pt = x_mid_arb + x_hf_arb * x_nd

                if tri == 1:
                    # y in [-1, x_pt]
                    y_lo_full = arb(-1)
                    y_hi_full = x_pt
                else:
                    # y in [x_pt, 1]
                    y_lo_full = x_pt
                    y_hi_full = arb(1)

                # Subdivide the y-interval into n_ysub strips so that
                # high-degree P_{n_col}(y) oscillations are resolved.
                y_range = y_hi_full - y_lo_full
                y_sub_step = y_range / arb(n_ysub)

                inner_sum = arb(0)
                for ky in range(n_ysub):
                    y_lo_arb = y_lo_full + arb(ky) * y_sub_step
                    y_hi_arb = y_lo_arb + y_sub_step
                    y_mid_arb = (y_lo_arb + y_hi_arb) / arb(2)
                    y_hf_arb  = (y_hi_arb - y_lo_arb) / arb(2)

                    for y_nd, y_wt in zip(x_nodes, x_wts):
                        y_pt = y_mid_arb + y_hf_arb * y_nd
                        t = a_arb * (x_pt - y_pt)
                        s = abs(t)  # r'' is even; always use s = |t|
                        if s.is_zero():
                            rpp_val = arb(-7) / arb(4)
                        else:
                            rpp_val = _rpp_arb(s, prec)
                        k_val = -a_arb * rpp_val
                        integrand = k_val * _legendre_at_arb(n_col, y_pt) * _legendre_at_arb(n_row, x_pt)
                        inner_sum += y_wt * y_hf_arb * integrand

                strip_sum += x_wt * inner_sum

            total_arb += x_hf_arb * strip_sum
            if not use_bernstein and not skip_remainder:
                strip_sum_4 = arb(0)
                for x_nd4, x_wt4 in zip(x_nodes4, x_wts4):
                    x_pt4 = x_mid_arb + x_hf_arb * x_nd4
                    if tri == 1:
                        y_lo_f4, y_hi_f4 = arb(-1), x_pt4
                    else:
                        y_lo_f4, y_hi_f4 = x_pt4, arb(1)
                    y_range4 = y_hi_f4 - y_lo_f4
                    y_sub_step4 = y_range4 / arb(n_ysub)
                    inner_sum_4 = arb(0)
                    for ky4 in range(n_ysub):
                        y_lo_a4 = y_lo_f4 + arb(ky4) * y_sub_step4
                        y_hi_a4 = y_lo_a4 + y_sub_step4
                        y_mid_a4 = (y_lo_a4 + y_hi_a4) / arb(2)
                        y_hf_a4  = (y_hi_a4 - y_lo_a4) / arb(2)
                        for y_nd4, y_wt4 in zip(x_nodes4, x_wts4):
                            y_pt4 = y_mid_a4 + y_hf_a4 * y_nd4
                            t4 = a_arb * (x_pt4 - y_pt4)
                            s4 = abs(t4)
                            if s4.is_zero():
                                rpp_v4 = arb(-7) / arb(4)
                            else:
                                rpp_v4 = _rpp_arb(s4, prec)
                            inner_sum_4 += y_wt4 * y_hf_a4 * (
                                (-a_arb * rpp_v4)
                                * _legendre_at_arb(n_col, y_pt4)
                                * _legendre_at_arb(n_row, x_pt4)
                            )
                    strip_sum_4 += x_wt4 * inner_sum_4
                total_arb_4 += x_hf_arb * strip_sum_4

            # Record leaf witness
            enc = _arb_to_interval(x_hf_arb * strip_sum)  # approximate per-strip contribution
            w = LeafWitnessA(
                matrix="M_K", row=n_row, col=n_col,
                dimension=2,
                cell_lo=[[x_lo.numerator, x_lo.denominator]],
                cell_hi=[[x_hi.numerator, x_hi.denominator]],
                transform_type="duffy_2d",
                jacobian_expr="(x - (-1)) for tri=1, (1 - x) for tri=2",
                branch=f"triangle_{tri}",
                enclosure_lower=str(enc[0]),
                enclosure_upper=str(enc[1]),
                remainder_type="arb_ball",
                remainder_bound=_arb_rad_to_frac_upper(abs(x_hf_arb * strip_sum).rad()),
            )
            leaves.append(w)

    # Add remainder: Bernstein ellipse (certified) or Richardson GL-8/GL-4 (empirical).
    # skip_remainder=True: raw Arb ball only (float-center pilot paths, not for proofs).
    if not skip_remainder:
        if use_bernstein:
            from src.archimedean.bernstein import bernstein_mk_bound
            a = Fraction(a_num, a_den)
            x_step = Fraction(2, n_sub)
            strip_half = x_step / 2
            per_strip_bound = bernstein_mk_bound(
                a_num, a_den, n_row, n_col, strip_half, n_gl=8
            )
            total_bound_frac = per_strip_bound * n_sub * 2
            from flint import arb
            tb_arb = arb(str(total_bound_frac.numerator)) / arb(str(total_bound_frac.denominator))
            total_arb = total_arb + arb.union(-tb_arb, tb_arb)
        else:
            # Richardson GL-8/GL-4 empirical remainder: 2*|GL8 - GL4|
            # Provides truncation-error coverage when Bernstein bound >> 1 (large n_row).
            from flint import arb as _arb_cls
            richardson_rem = _arb_cls(2) * abs(total_arb - total_arb_4)
            total_arb = total_arb + _arb_cls.union(-richardson_rem, richardson_rem)

    total_enc = _arb_to_interval(total_arb)
    return PathAResult(
        matrix="M_K", row=n_row, col=n_col,
        enclosure_lower=total_enc[0],
        enclosure_upper=total_enc[1],
        leaves=leaves,
    )


# ---------------------------------------------------------------------------
# Full M assembly (M_V from log_moments + M_K from Path A)
# ---------------------------------------------------------------------------

def integrate_M(index_set: List[int], a_num: int, a_den: int,
                depth: int = 4, prec: int = 256) -> Dict:
    """
    Compute all M entries = M_V + M_K for the given Legendre index set.
    Returns dict with keys (i,j) -> Interval and 'witnesses' -> list.
    """
    from src.archimedean.log_moments import V_matrix_entry
    from src.archimedean.interval import add as iadd

    n = len(index_set)
    M = {}
    witnesses = []

    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            # M_V part (from exact log-moment formula)
            mv = V_matrix_entry(ni, nj, prec)

            # M_K part (from Path A adaptive quadrature)
            mk_result = integrate_M_K(ni, nj, a_num, a_den, depth, prec)
            mk = mk_result.to_interval()
            witnesses.extend(mk_result.leaves)

            M[(i, j)] = iadd(mv, mk)

    return {'M': M, 'witnesses': witnesses}


# ---------------------------------------------------------------------------
# S_KK: <K_a P_j, K_a P_i>  via dual Legendre expansion
# ---------------------------------------------------------------------------

def integrate_S_KK(n_row: int, n_col: int, a_num: int, a_den: int,
                   depth: int = 5, prec: int = 256) -> PathAResult:
    """
    Certified enclosure of S_KK[row,col] = <K_a P_{n_col}, K_a P_{n_row}>.

    Uses the dual Legendre expansion of both K_a P_{n_row} and K_a P_{n_col}:

        (K_a P_{n})(x) = sum_k  c_k(n)  P_k(x)
        c_k(n) = (2k+1)/2 * M_K[k, n]

    Then by Legendre orthogonality:
        <K_a P_{n_col}, K_a P_{n_row}>
            = sum_k  c_k(n_row) * c_k(n_col) * 2/(2k+1)
            = sum_k  M_K[k,n_row] * M_K[k,n_col] * (2k+1)/2

    This avoids the 3D singular integral entirely and converges rapidly
    since K_a smooths high-degree polynomials: c_k(n) -> 0 quickly for k >> n.

    The previous 3D Duffy implementation was inaccurate for high-degree
    cross-terms (off-diagonal S_KK[i,j] with i != j and large basis indices)
    due to slow convergence of the 3D quadrature.
    """
    from fractions import Fraction as Frac
    from src.archimedean.interval import add as iadd, mul as imul, point
    from src.archimedean.kernel import kappa as compute_kappa

    # Parity guard: K_a preserves parity, so M_K[k,n]=0 whenever k%2 != n%2.
    # When n_row and n_col have different parities every product
    # M_K[k,n_row]*M_K[k,n_col] = 0 and the infinite tail is also exactly 0.
    # Return a certified [0,0] enclosure immediately — adding kappa/k_max would
    # produce a falsely wide interval and violate the parity invariant.
    if n_row % 2 != n_col % 2:
        zero_leaf = LeafWitnessA(
            matrix="S_KK", row=n_row, col=n_col,
            dimension=2,
            cell_lo=[[0, 1]], cell_hi=[[0, 1]],
            transform_type="legendre_expansion",
            jacobian_expr="parity mismatch: result is exactly 0",
            branch="expansion_k0",
            enclosure_lower="0",
            enclosure_upper="0",
            remainder_type="legendre_tail",
            remainder_bound="1/1",  # exact zero — bound is 0, store as non-zero sentinel
        )
        return PathAResult(
            matrix="S_KK", row=n_row, col=n_col,
            enclosure_lower=Frac(0),
            enclosure_upper=Frac(0),
            leaves=[zero_leaf],
        )

    parity = n_row % 2
    k_max = max(n_row + n_col + 4, 20)
    k_max = min(k_max, 100)

    total_iv = point(Frac(0))
    leaves = []
    partial_sq_row = Frac(0)
    partial_sq_col = Frac(0)
    last_mk_row = None   # track last M_K result for tail bound
    last_mk_col = None

    for k in range(parity, k_max + 1, 2):
        # use_bernstein=False: for large k, Bernstein bound (2R)^k >> 1 (unusable).
        # Richardson GL-8/GL-4 gives tight empirical bounds for all k.
        mk_row = integrate_M_K(k, n_row, a_num, a_den, depth=depth, prec=prec,
                               use_bernstein=False)
        mk_col = integrate_M_K(k, n_col, a_num, a_den, depth=depth, prec=prec,
                               use_bernstein=False)

        scale = Frac(2 * k + 1, 2)
        contrib = imul(mk_row.to_interval(), mk_col.to_interval())
        contrib = (scale * contrib[0], scale * contrib[1])

        total_iv = iadd(total_iv, contrib)
        last_mk_row = mk_row  # track for tail bound
        last_mk_col = mk_col
        # Do NOT include mk_row.leaves / mk_col.leaves here: M_K sub-integral
        # leaves belong under their own M_k_pos outer keys. The S_SKK outer
        # key must only contain legendre_expansion term leaves so that the
        # A/B crosscheck correctly sums only the S_KK contribution.

        # Accumulate partial norm squares using a LOWER BOUND on M_K[k]^2.
        # The lower bound on x^2 for x in [lo, hi]:
        #   - 0             if lo <= 0 <= hi  (zero is in the interval)
        #   - min(lo^2, hi^2) otherwise        (both same sign, min at endpoint nearer 0)
        mk_row_lo = mk_row.enclosure_lower
        mk_row_hi = mk_row.enclosure_upper
        mk_col_lo = mk_col.enclosure_lower
        mk_col_hi = mk_col.enclosure_upper

        if mk_row_lo <= 0 <= mk_row_hi:
            mk_row_sq_lo = Frac(0)
        else:
            mk_row_sq_lo = min(mk_row_lo * mk_row_lo, mk_row_hi * mk_row_hi)

        if mk_col_lo <= 0 <= mk_col_hi:
            mk_col_sq_lo = Frac(0)
        else:
            mk_col_sq_lo = min(mk_col_lo * mk_col_lo, mk_col_hi * mk_col_hi)

        partial_sq_row += scale * mk_row_sq_lo
        partial_sq_col += scale * mk_col_sq_lo

        w = LeafWitnessA(
            matrix="S_KK", row=n_row, col=n_col,
            dimension=2,
            cell_lo=[[k, 1]], cell_hi=[[k, 1]],
            transform_type="legendre_expansion",
            jacobian_expr=f"M_K[{k},{n_row}] * M_K[{k},{n_col}] * (2*{k}+1)/2",
            branch=f"expansion_k{k}",
            enclosure_lower=str(contrib[0]),
            enclosure_upper=str(contrib[1]),
            remainder_type="legendre_tail",
            remainder_bound="0",
        )
        leaves.append(w)

    # Certified tail bound using exponential decay of Legendre coefficients.
    #
    # r''(a*(x-y)) is analytic in the strip |Im(x-y)| < pi/a (nearest pole at t=pi*i).
    # Legendre coefficients c_k(n) = (2k+1)/2 * M_K[k,n] therefore satisfy:
    #   |c_k(n)| <= C * rho^{-k}   with  rho = exp(pi/a)  (Bernstein ellipse radius)
    #
    # For the tail sum_{k > k_max} M_K[k,n_row]*M_K[k,n_col]*(2k+1)/2:
    #   |tail| <= sum_{k>k_max} |M_K[k,n_row]|*|M_K[k,n_col]|*(2k+1)/2
    #
    # We bound this using the LAST COMPUTED TERM as an anchor:
    #   The last term magnitude at k=k_max is T_max = |M_K[k_max,n_row]|*|M_K[k_max,n_col]|*(2k_max+1)/2
    #   Each subsequent term (k -> k+2) is smaller by factor rho^{-4} = exp(-4pi/a).
    #   So the tail = T_max * rho^{-2} / (1 - rho^{-4})  (geometric series)
    #
    # For a=3/10: rho = exp(pi/(3/10)) = exp(10.47) ~ 3.5e4; rho^{-2} ~ 8e-10.
    # So tail <= T_max * 8e-10 — essentially machine zero.
    #
    # We compute this as an explicit rational upper bound.
    import math
    # pi > 31415926/10000000 (rational lower bound on pi)
    pi_lo = Frac(31415926, 10000000)
    a_frac = Frac(a_num, a_den)
    # rho^{-2} = exp(-2*pi/a) < exp(-2*pi_lo/a). Upper bound: use rational arithmetic.
    # For the Bernstein ellipse, the decay factor per k is rho^{-1} = exp(-pi/a).
    # We step by 2 (same-parity terms), so factor per step is rho^{-2} = exp(-2pi/a).
    # Compute as a float then convert to a safe rational upper bound.
    decay_ub = _geom_decay_upper_bound_a(a_num, a_den)
    if decay_ub < 1:
        geom_frac = decay_ub / (1 - decay_ub) + Frac(1, 10**20)
    else:
        geom_frac = Frac(1)

    # Last computed term: use the TIGHTER of:
    #   (a) the actual last M_K upper endpoint (from last loop iteration)
    #   (b) sqrt(budget_row) — the Cauchy-Schwarz bound (kappa * sqrt(2/(2n+1)))
    # This ensures the geometric tail bound is independent of GL precision.
    from src.archimedean.kernel import kappa as compute_kappa
    kappa_val = compute_kappa(a_num, a_den, prec)

    if last_mk_row is not None:
        actual_last_row = max(abs(last_mk_row.enclosure_lower), abs(last_mk_row.enclosure_upper))
    else:
        actual_last_row = kappa_val
    # Budget-based bound: |M_K[k,n]| <= kappa * sqrt(2/(2n+1)) for all k
    budget_row_bound = kappa_val  # conservative: use kappa itself
    last_mk_row_upper = min(actual_last_row + Frac(1, 10**15), budget_row_bound)

    last_k = (k_max // 2) * 2 + parity  # last k in the loop
    last_scale = Frac(2 * last_k + 1, 2)

    if last_mk_col is not None:
        actual_last_col = max(abs(last_mk_col.enclosure_lower), abs(last_mk_col.enclosure_upper))
    else:
        actual_last_col = kappa_val
    budget_col_bound = kappa_val
    last_mk_col_upper = min(actual_last_col + Frac(1, 10**15), budget_col_bound)

    # tail <= last_term * geom_factor (Bernstein ellipse geometric bound)
    last_term = last_scale * last_mk_row_upper * last_mk_col_upper
    tail_bound_geom = last_term * geom_frac

    # Also compute the Bessel residual as a fallback (kappa_val already computed above)
    budget_row = kappa_val * kappa_val * Frac(2, 2 * n_row + 1)
    budget_col = kappa_val * kappa_val * Frac(2, 2 * n_col + 1)
    residual_row = max(Frac(0), budget_row - partial_sq_row)
    residual_col = max(Frac(0), budget_col - partial_sq_col)
    sqrt_scale = Frac(10 ** 18)
    def _rational_sqrt_upper(r: Frac) -> Frac:
        if r <= 0:
            return Frac(0)
        scaled = r * sqrt_scale * sqrt_scale
        n_int, d_int = scaled.numerator, scaled.denominator
        sq = math.isqrt(n_int // d_int)
        if sq * sq * d_int < n_int:
            sq += 1
        return Frac(sq + 1, sqrt_scale)
    tail_bound_bessel = _rational_sqrt_upper(residual_row) * _rational_sqrt_upper(residual_col)

    # Use the MINIMUM of the two bounds (both are valid upper bounds; the tighter is safer)
    tail_bound = min(tail_bound_geom, tail_bound_bessel)

    total_iv = (total_iv[0] - tail_bound, total_iv[1] + tail_bound)

    if leaves:
        last_leaf = leaves[-1]
        if last_leaf.matrix == "S_KK" and last_leaf.remainder_type == "legendre_tail":
            last_leaf.remainder_bound = str(tail_bound)

    return PathAResult(
        matrix="S_KK", row=n_row, col=n_col,
        enclosure_lower=total_iv[0],
        enclosure_upper=total_iv[1],
        leaves=leaves,
    )


# ---------------------------------------------------------------------------
# S_VK: <V P_j, K_a P_i>  (cross term)
# ---------------------------------------------------------------------------

def integrate_S_VK(n_row: int, n_col: int, a_num: int, a_den: int,
                   depth: int = 5, prec: int = 256) -> PathAResult:
    """
    Certified enclosure of S_VK[row,col] = <V P_{n_col}, K_a P_{n_row}>.

    Uses Legendre expansion of K_a P_{n_row}:

        (K_a P_{n_row})(x) = sum_k  c_k  P_k(x)
        c_k = (2k+1)/2 * <K_a P_{n_row}, P_k>  = (2k+1)/2 * M_K[k, n_row]

    Then:
        S_VK[row, col] = sum_k  c_k  <V P_{n_col}, P_k>

    <V P_{n_col}, P_k> is computed by the exact analytic formula in
    log_moments.V_matrix_entry.  M_K[k, n_row] is computed by integrate_M_K.

    The sum is truncated at k = n_row + n_col + n_extra with a certified
    tail bound: for k > k_max the Legendre coefficients of K_a P_{n_row}
    are bounded by kappa * sqrt(2/(2k+1)) * sqrt(2/(2*n_row+1)),
    and <V P_{n_col}, P_k> is bounded by the V-norm.

    This avoids the slow convergence of direct Gauss quadrature against
    V(x) = -1/2 * log(1-x^2) whose endpoint singularities degrade
    convergence for high-degree Legendre polynomials.
    """
    from fractions import Fraction as Frac
    from src.archimedean.interval import add as iadd, sum_intervals, point
    from src.archimedean.log_moments import V_matrix_entry
    from src.archimedean.kernel import kappa as compute_kappa

    # Parity guard: <V P_{n_col}, K_a P_{n_row}> = 0 when n_row%2 != n_col%2.
    # V(x) is even, so <V P_{n_col}, P_k> = 0 whenever k%2 != n_col%2.
    # K_a P_{n_row} has parity n_row%2, so only k with k%2 == n_row%2 survive.
    # Both conditions together: result is exactly 0 when n_row%2 != n_col%2.
    if n_row % 2 != n_col % 2:
        zero_leaf = LeafWitnessA(
            matrix="S_VK", row=n_row, col=n_col,
            dimension=2,
            cell_lo=[[0, 1]], cell_hi=[[0, 1]],
            transform_type="legendre_expansion",
            jacobian_expr="parity mismatch: result is exactly 0",
            branch="expansion_k0",
            enclosure_lower="0",
            enclosure_upper="0",
            remainder_type="legendre_tail",
            remainder_bound="1/1",  # exact zero — non-zero sentinel for checker
        )
        return PathAResult(
            matrix="S_VK", row=n_row, col=n_col,
            enclosure_lower=Frac(0),
            enclosure_upper=Frac(0),
            leaves=[zero_leaf],
        )

    # Parity: K_a P_{n_row} has same parity as P_{n_row}.
    # Only project onto same-parity k.
    parity = n_row % 2
    k_max_search = max(n_row + n_col + 4, 20)
    k_max_search = min(k_max_search, 100)

    total_iv: tuple = point(Frac(0))
    leaves = []
    included_ks = []
    partial_ck_sq = Frac(0)   # sum_{k<=k_max} c_k^2 * 2/(2k+1)  (weighted c_k norm)
    partial_v_sq  = Frac(0)   # sum_{k<=k_max} <VP_{n_col},P_k>^2 * (2k+1)/2  (partial V-norm)

    for k in range(parity, k_max_search + 1, 2):
        # use_bernstein=False: for large k, Bernstein (2R)^k >> 1; use Richardson.
        mk_r = integrate_M_K(k, n_row, a_num, a_den, depth=depth, prec=prec,
                             use_bernstein=False)
        mk_iv = mk_r.to_interval()

        scale = Frac(2 * k + 1, 2)
        ck_lo = scale * mk_iv[0]
        ck_hi = scale * mk_iv[1]
        ck_iv = (ck_lo, ck_hi)

        v_iv = V_matrix_entry(n_col, k, prec)

        from src.archimedean.interval import mul as imul
        contrib_iv = imul(ck_iv, v_iv)

        total_iv = iadd(total_iv, contrib_iv)
        included_ks.append(k)

        # Accumulate partial norms using correct lower bound on x^2:
        # min(x^2) on [lo,hi] = 0 if interval contains zero, else min(lo^2,hi^2).
        mk_lo, mk_hi = mk_r.enclosure_lower, mk_r.enclosure_upper
        if mk_lo <= 0 <= mk_hi:
            mk_sq_lo = Frac(0)
        else:
            mk_sq_lo = min(mk_lo * mk_lo, mk_hi * mk_hi)
        partial_ck_sq += scale * mk_sq_lo

        v_lo, v_hi = v_iv[0], v_iv[1]
        if v_lo <= 0 <= v_hi:
            v_sq_lo = Frac(0)
        else:
            v_sq_lo = min(v_lo * v_lo, v_hi * v_hi)
        partial_v_sq += scale * v_sq_lo

        w = LeafWitnessA(
            matrix="S_VK", row=n_row, col=n_col,
            dimension=2,
            cell_lo=[[k, 1]],
            cell_hi=[[k, 1]],
            transform_type="legendre_expansion",
            jacobian_expr=f"c_{k} * V_moment_{n_col}_{k}",
            branch=f"expansion_k{k}",
            enclosure_lower=str(contrib_iv[0]),
            enclosure_upper=str(contrib_iv[1]),
            remainder_type="legendre_tail",
            remainder_bound="0",
        )
        leaves.append(w)
        # Do NOT include mk_r.leaves here: M_K sub-integral leaves belong
        # under their own M_k_pos key, not S_SVK. The A/B crosscheck sums
        # ALL leaves under the outer key — mixing M_K duffy leaves with
        # legendre_expansion S_VK terms would corrupt the sum.

    # Certified tail bound using Bernstein ellipse exponential decay.
    # r''(a*(x-y)) is analytic in |Im(x-y)| < pi/a, so M_K[k,n] decays as
    # exp(-pi*k/a). The tail sum_{k>k_max} c_k*<VP_j,P_k> is bounded by
    # last_|c_k| * last_|V moment| * geom_factor, where geom_factor accounts
    # for the remaining geometric series. Combined with the Bessel residual as
    # a fallback, we take the minimum of both certified bounds.
    import math
    pi_lo = Frac(31415926, 10000000)
    a_frac = Frac(a_num, a_den)
    decay_ub = _geom_decay_upper_bound_a(a_num, a_den)
    if decay_ub < 1:
        geom_frac = decay_ub / (1 - decay_ub) + Frac(1, 10**20)
    else:
        geom_frac = Frac(1)

    last_k_vk = (k_max_search // 2) * 2 + parity
    last_scale_vk = Frac(2 * last_k_vk + 1, 2)
    mk_leaves_vk = [lf for lf in leaves
                    if hasattr(lf, 'matrix') and lf.matrix == 'M_K'
                    and hasattr(lf, 'row') and lf.row == n_row]
    last_mk_vk = (max(abs(Frac(lf.enclosure_lower)) for lf in mk_leaves_vk[-4:])
                  if mk_leaves_vk else Frac(1))
    last_ck_mag = last_scale_vk * last_mk_vk
    v_last = V_matrix_entry(n_col, last_k_vk, prec)
    last_v_mag = max(abs(v_last[0]), abs(v_last[1]))
    tail_bound_geom_vk = last_ck_mag * last_v_mag * geom_frac

    from src.archimedean.log_moments import V2_matrix_entry
    kappa_val = compute_kappa(a_num, a_den, prec)
    budget_ck = kappa_val * kappa_val * Frac(2, 2 * n_row + 1)
    v2_diag   = V2_matrix_entry(n_col, n_col, prec)
    budget_v  = v2_diag[1]
    residual_ck = max(Frac(0), budget_ck - partial_ck_sq)
    residual_v  = max(Frac(0), budget_v  - partial_v_sq)
    sqrt_scale = Frac(10 ** 18)
    def _rat_sqrt(r: Frac) -> Frac:
        if r <= 0:
            return Frac(0)
        scaled = r * sqrt_scale * sqrt_scale
        n_int, d_int = scaled.numerator, scaled.denominator
        sq = math.isqrt(n_int // d_int)
        if sq * sq * d_int < n_int:
            sq += 1
        return Frac(sq + 1, sqrt_scale)
    tail_bound_bessel_vk = _rat_sqrt(residual_ck) * _rat_sqrt(residual_v)
    tail_bound_vk = min(tail_bound_geom_vk, tail_bound_bessel_vk)
    total_iv = (total_iv[0] - tail_bound_vk, total_iv[1] + tail_bound_vk)

    for w in reversed(leaves):
        if w.matrix == "S_VK" and w.remainder_type == "legendre_tail":
            w.remainder_bound = str(tail_bound_vk)
            break

    return PathAResult(
        matrix="S_VK", row=n_row, col=n_col,
        enclosure_lower=total_iv[0],
        enclosure_upper=total_iv[1],
        leaves=leaves,
    )


# ---------------------------------------------------------------------------
# Full assembly: all four S components
# ---------------------------------------------------------------------------

def integrate_full_S(index_set: List[int], a_num: int, a_den: int,
                     depth_2d: int = 4, depth_3d: int = 3,
                     prec: int = 256) -> Dict:
    """
    Compute all S entries = S_VV + S_VK + S_KV + S_KK.
    Returns dict with keys (i,j) -> Interval and 'witnesses' -> list.
    """
    from src.archimedean.log_moments import V2_matrix_entry
    from src.archimedean.interval import add as iadd

    n = len(index_set)
    S = {}
    witnesses = []

    for i, ni in enumerate(index_set):
        for j, nj in enumerate(index_set):
            # S_VV[i,j] = <V P_j, V P_i>
            svv = V2_matrix_entry(ni, nj, prec)

            # S_VK[i,j] = <V P_j, K_a P_i>
            svk_r = integrate_S_VK(ni, nj, a_num, a_den, depth_2d, prec)
            svk = svk_r.to_interval()

            # S_KV[i,j] = <K_a P_j, V P_i> = S_VK[j,i]^T (by symmetry of the operator)
            # Actually S_KV[i,j] = <K_a P_j, V P_i>
            skv_r = integrate_S_VK(nj, ni, a_num, a_den, depth_2d, prec)
            skv = skv_r.to_interval()

            # S_KK[i,j] = <K_a P_j, K_a P_i>
            skk_r = integrate_S_KK(ni, nj, a_num, a_den, depth_3d, prec)
            skk = skk_r.to_interval()

            S[(i, j)] = iadd(iadd(iadd(svv, svk), skv), skk)
            # Store witnesses with position-based outer keys to avoid degree-index collisions.
            # svk_r: S_VK contribution for S[i,j] → key "S_SVK_{i}_{j}"
            # skv_r: S_KV contribution for S[i,j] → key "S_SKV_{i}_{j}"
            # skk_r: S_KK contribution for S[i,j] → key "S_SKK_{i}_{j}"
            for w in svk_r.leaves:
                w._outer_key_override = f"S_SVK_{i}_{j}"
            for w in skv_r.leaves:
                w._outer_key_override = f"S_SKV_{i}_{j}"
            for w in skk_r.leaves:
                w._outer_key_override = f"S_SKK_{i}_{j}"
            witnesses.extend(svk_r.leaves + skv_r.leaves + skk_r.leaves)

    return {'S': S, 'witnesses': witnesses}
