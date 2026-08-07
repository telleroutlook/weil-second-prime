"""
Endpoint log-moment integral library.

Generates certified interval enclosures for all matrix entries involving the
potential V(x) = -1/2 * log(1 - x^2) and its square V^2(x).

The generating family:
    I_m(alpha, beta) = integral_{-1}^{1} x^m (1-x)^alpha (1+x)^beta dx

Differentiate w.r.t. alpha and/or beta at (0,0) to obtain:
    dI/dalpha  -> integral of x^m * log(1-x)
    dI/dbeta   -> integral of x^m * log(1+x)
    d^2I/dalpha^2 -> integral of x^m * log^2(1-x)
    d^2I/dalphadbeta -> integral of x^m * log(1-x)*log(1+x)

The Beta function identity:
    I_m(alpha, beta) = sum_{k=0}^{m} C(m,k) * (-1)^k * B(k+alpha+1, beta+1)
with B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b).

For integer a,b: B(a,b) = (a-1)!(b-1)! / (a+b-1)!
The derivative dB/dalpha at alpha=0 involves the digamma function psi(a) = H_{a-1} - gamma.

Reference: 28-day plan Section 6.1; proof/legendre-tail.tex Section 2.
D5 deliverable: log-moments.json
"""

from fractions import Fraction
from typing import List, Tuple

Interval = Tuple[Fraction, Fraction]


# ---------------------------------------------------------------------------
# Exact Beta function and its derivatives at integer arguments
# ---------------------------------------------------------------------------

def _factorial(n: int) -> int:
    result = 1
    for k in range(2, n + 1):
        result *= k
    return result


def beta_exact(a: int, b: int) -> Fraction:
    """
    B(a, b) = (a-1)! * (b-1)! / (a+b-1)!  for positive integers a, b.
    Returns exact Fraction.
    """
    return Fraction(_factorial(a - 1) * _factorial(b - 1), _factorial(a + b - 1))


def Im_exact(m: int, alpha: int, beta: int) -> Fraction:
    """
    I_m(alpha, beta) = integral_{-1}^{1} x^m (1-x)^alpha (1+x)^beta dx
    for non-negative integers alpha, beta.

    Uses: x^m (1-x)^alpha = sum_{k=0}^{m} C(m,k)(-1)^k x^k (1-x)^alpha ...
    Actually use the substitution x = 2t-1 to map to [0,1]:
        I_m(alpha,beta) = 2^{m+alpha+beta+1} * B(m+alpha+beta+...  )

    Correct direct formula via binomial on (1-x)^alpha = sum_j C(a,j)(-x)^j:
        I_m(alpha,beta) = sum_{j=0}^{alpha} C(alpha,j)(-1)^j * J(m+j, beta)
    where J(p, beta) = integral_{-1}^{1} x^p (1+x)^beta dx
                     = sum_{k=0}^{beta} C(beta,k) * integral_{-1}^{1} x^{p+k} dx
                     = sum_{k=0}^{beta} C(beta,k) * [2/(p+k+1) if p+k even, 0 otherwise]
    """
    from math import comb
    total = Fraction(0)
    for j in range(alpha + 1):
        cj = Fraction(comb(alpha, j) * ((-1) ** j))
        p = m + j
        for k in range(beta + 1):
            ck = Fraction(comb(beta, k))
            pk = p + k
            if pk % 2 == 0:  # integral of x^pk from -1 to 1 = 2/(pk+1)
                total += cj * ck * Fraction(2, pk + 1)
    return total


# ---------------------------------------------------------------------------
# Derivatives using Arb (for log terms)
# ---------------------------------------------------------------------------

def _arb_to_interval(x) -> Interval:
    """Convert arb ball to outward-rounded Fraction interval."""
    digits = 60
    M, R, E = x.mid_rad_10exp(digits)
    M, R, E = int(M), int(R), int(E)
    if M == 0 and R == 0:
        return Fraction(0), Fraction(0)
    if E >= 0:
        scale = Fraction(10 ** E)
    else:
        scale = Fraction(1, 10 ** (-E))
    mid = Fraction(M) * scale
    rad = Fraction(R) * scale
    ulp = abs(scale)
    return mid - rad - ulp, mid + rad + ulp


def _legendre_coeffs(n: int) -> List[Fraction]:
    """
    Return coefficients [c_0, c_1, ..., c_n] such that
    P_n(x) = sum_k c_k * x^k, using the standard Legendre recurrence.
    """
    if n == 0:
        return [Fraction(1)]
    if n == 1:
        return [Fraction(0), Fraction(1)]
    p_prev = [Fraction(1)]
    p_curr = [Fraction(0), Fraction(1)]
    for k in range(2, n + 1):
        # P_k = ((2k-1)*x*P_{k-1} - (k-1)*P_{k-2}) / k
        # x * p_curr: shift coefficients up by one
        x_p_curr = [Fraction(0)] + p_curr
        p_next = [Fraction(0)] * (k + 1)
        for i in range(k + 1):
            xc = x_p_curr[i] if i < len(x_p_curr) else Fraction(0)
            pc = p_prev[i] if i < len(p_prev) else Fraction(0)
            p_next[i] = (Fraction(2 * k - 1) * xc - Fraction(k - 1) * pc) / k
        p_prev = p_curr
        p_curr = p_next
    return p_curr


def V_matrix_entry(n_row: int, n_col: int, prec: int = 256) -> Interval:
    """
    <V P_{n_col}, P_{n_row}> where V(x) = -1/2 * log(1-x^2).

    = -1/2 * integral_{-1}^{1} log(1-x^2) * P_{n_row}(x) * P_{n_col}(x) dx
    = -1/2 * (integral log(1-x) * P_n_row * P_n_col dx
              + integral log(1+x) * P_n_row * P_n_col dx)

    Computed via Arb numerical integration at high precision.
    """
    from flint import arb, ctx
    ctx.prec = prec

    # Expand P_{n_row}(x) * P_{n_col}(x) in monomial basis
    c_row = _legendre_coeffs(n_row)
    c_col = _legendre_coeffs(n_col)

    # Product polynomial coefficients
    deg = n_row + n_col
    prod = [Fraction(0)] * (deg + 1)
    for i, ci in enumerate(c_row):
        for j, cj in enumerate(c_col):
            prod[i + j] += ci * cj

    # Compute integral_{-1}^{1} x^m * log(1-x^2) dx for each monomial
    # = integral log(1-x)*x^m dx + integral log(1+x)*x^m dx
    # Odd powers: both integrals cancel by symmetry -> 0 for odd m
    # Even powers: both equal by symmetry -> 2 * integral_0^1 x^m log(1-x^2) dx

    # Use Arb for the log integrals
    total = arb(0)
    for m, cm in enumerate(prod):
        if cm == 0:
            continue
        # integral_{-1}^{1} x^m * log(1 - x^2) dx
        # = 2 * integral_0^1 x^m * log(1-x^2) dx   (if m even, else 0)
        if m % 2 != 0:
            continue
        cm_arb = arb(str(cm.numerator)) / arb(str(cm.denominator))
        # Compute via Beta derivatives:
        # integral_0^1 x^m * log(1-x^2) dx = d/dalpha [B((m+1)/2, alpha+1/2)] / 2 at alpha=0
        # Simpler: use Arb quad or the known formula
        # integral_0^1 x^m * log(1-x) dx = -(H_m - H_0 + ...) via digamma
        # We use: d/da B(a,b)|_{a=integer} involves polygamma
        # Numerically:
        # integral_{-1}^{1} x^m * log(1-x^2) dx at m=0:
        #   = integral_{-1}^1 log(1-x^2) dx = -4 + 2*log(4) ... known
        # Use Arb integration
        from flint import arb
        # arb.integral is not available in python-flint 0.9; use quadrature manually
        # Build the integrand symbolically via polygamma
        # integral_0^1 t^m log(1-t^2) dt: substitute t^2 = u -> 1/2 integral_0^1 u^{(m-1)/2} log(1-u) du
        # = 1/2 * d/db [B((m+1)/2, b+1)] at b=0
        # B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)
        # d/db log B = psi(b) - psi(a+b)  ->  at b=1: psi(1)-psi(a+1) = -gamma - (H_a - gamma) = -H_a
        # Wait: d/db B(a,b)|_{b=1} = B(a,1) * (psi(1) - psi(a+1)) = (1/a)*(-H_a)  [since psi(n+1)=H_n-gamma, psi(1)=-gamma]
        # So integral_0^1 t^m log(1-t^2) dt = 1/2 * (1/((m+1)/2)) * (-H_{(m+1)/2}) ...
        # This only works when (m+1)/2 is integer, i.e., m odd — but we only have even m here.
        # For even m: let m = 2p. integral_0^1 t^{2p} log(1-t^2) dt
        # sub t^2=u: (1/2) integral_0^1 u^{p-1/2} log(1-u) du = (1/2) d/db B(p+1/2, b+1)|_{b=0}
        # This involves half-integer Gamma — use Arb directly.
        half_m_plus_1 = arb(m + 1) / arb(2)   # (m+1)/2, half-integer for even m
        # integral_0^1 u^{(m-1)/2} log(1-u) du = d/db B((m+1)/2, b+1) at b=0
        # = B((m+1)/2, 1) * (psi(1) - psi((m+1)/2 + 1))
        # B(a,1) = 1/a; psi(1) = -gamma; psi(a+1) = psi(a) + 1/a
        # For half-integer a = n+1/2: psi(n+1/2) = -gamma - 2*log(2) + 2*(1 + 1/3 + ... + 1/(2n-1))
        # Use Arb polygamma directly
        psi1  = -arb.const_euler()   # psi(1) = -gamma
        psi_a1 = arb.digamma(half_m_plus_1 + arb(1))  # psi((m+1)/2 + 1)
        B_a1  = arb(1) / half_m_plus_1  # B((m+1)/2, 1) = 1/a
        inner = B_a1 * (psi1 - psi_a1)   # d/db B at b=0 (via chain rule * B * digamma diff)
        # integral_{-1}^1 x^m log(1-x^2) dx = 2 * integral_0^1 * (1/2) * inner = inner
        total += cm_arb * inner

    result_arb = -arb("1/2") * total
    return _arb_to_interval(result_arb)


def V2_matrix_entry(n_row: int, n_col: int, prec: int = 256) -> Interval:
    """
    <V P_{n_col}, V P_{n_row}> = integral_{-1}^1 V(x)^2 P_{n_row}(x) P_{n_col}(x) dx

    V(x)^2 = (1/4) * log^2(1-x^2)

    For even monomial x^m (m=2p), the integral reduces to:
        integral_{-1}^1 x^m log^2(1-x^2) dx
          = d^2/db^2 B((m+1)/2, b+1)|_{b=0}
          = B((m+1)/2, 1) * [(psi(1)-psi((m+1)/2+1))^2 + psi'(1) - psi'((m+1)/2+1)]

    psi(1) = -gamma, psi(n+1) = H_n - gamma  (exact harmonic number formula)
    psi'(1) = pi^2/6
    psi'(p+3/2) = pi^2/2 - 4*sum_{k=0}^p 1/(2k+1)^2  (analytic half-integer recurrence)

    NO finite-difference approximation is used. The digamma and polygamma at
    half-integer arguments are evaluated via these closed-form recurrences,
    with pi^2 evaluated via Arb at full precision.
    """
    from flint import arb, ctx
    ctx.prec = prec

    c_row = _legendre_coeffs(n_row)
    c_col = _legendre_coeffs(n_col)

    deg = n_row + n_col
    prod = [Fraction(0)] * (deg + 1)
    for i, ci in enumerate(c_row):
        for j, cj in enumerate(c_col):
            prod[i + j] += ci * cj

    total = arb(0)
    for m, cm in enumerate(prod):
        if cm == 0 or m % 2 != 0:
            continue
        cm_arb = arb(str(cm.numerator)) / arb(str(cm.denominator))
        # m = 2p, argument = (m+1)/2 = p + 1/2 is a half-integer.
        # half_m_plus_1 = (m+1)/2 = p + 1/2
        p = m // 2   # m = 2p, so (m+1)/2 = p + 1/2
        # B((m+1)/2, 1) = 1 / ((m+1)/2) = 2/(m+1)
        B_a1 = Fraction(2, m + 1)
        B_a1_arb = arb(str(B_a1.numerator)) / arb(str(B_a1.denominator))

        # psi(1) = -gamma  (exact via Arb)
        psi1 = -arb.const_euler()

        # psi((m+1)/2 + 1) = psi(p + 3/2)
        # Use recurrence: psi(p+3/2) = psi(1/2) + sum_{k=0}^{p} 1/(k+1/2)
        #                             = -gamma - 2*log(2) + 2*sum_{k=0}^{p} 1/(2k+1)
        # (standard half-integer digamma formula)
        h_sum = Fraction(0)
        for k in range(p + 1):
            h_sum += Fraction(1, 2 * k + 1)
        h_sum_arb = arb(str(h_sum.numerator)) / arb(str(h_sum.denominator))
        psi_a1 = -arb.const_euler() - 2 * arb.log(arb(2)) + 2 * h_sum_arb

        # psi'(1) = pi^2/6  (exact via Arb)
        psip1 = arb.pi() ** 2 / arb(6)

        # psi'(p+3/2) = pi^2/2 - 4*sum_{k=0}^{p} 1/(2k+1)^2  (analytic recurrence)
        sq_sum = Fraction(0)
        for k in range(p + 1):
            sq_sum += Fraction(1, (2 * k + 1) ** 2)
        sq_sum_arb = arb(str(sq_sum.numerator)) / arb(str(sq_sum.denominator))
        psip_a1 = arb.pi() ** 2 / arb(2) - arb(4) * sq_sum_arb

        diff_psi = psi1 - psi_a1
        d2_logB = diff_psi ** 2 + psip1 - psip_a1
        inner = B_a1_arb * d2_logB
        total += cm_arb * inner

    result_arb = arb("1/4") * total
    return _arb_to_interval(result_arb)
