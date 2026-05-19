"""
g0_j3.py — Golden octonionic order G_0 = I + I*ell and the hermitian
module J_3(G_0).

An octonion element is represented as a pair (a, b) with a, b in the
icosian ring I (Quat objects with Zphi coefficients).  The
Cayley-Dickson product is
    (a, b)(c, d) = (a c - conj(d) b, d a + b conj(c)),
and the involution is (a, b)^* = (conj(a), -b).

J_3(G_0) is the 27-dimensional O_K-module of hermitian 3x3 matrices
over G_0, with basis E_11, E_22, E_33 and F_ij(e_a) for the 24
off-diagonal slot/basis pairs.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple

from Z_phi import Zphi, ZERO, ONE, HALF, PHI, PHI_INV
from quaternion import Quat, ICOSIAN_BASIS, ICOSIAN_NAMES, QH_ONE


@dataclass(frozen=True)
class Oct:
    """Octonion (a, b) representing a + b * ell, a, b in H_K (Quat)."""
    a: Quat
    b: Quat

    def __add__(self, other):
        return Oct(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        return Oct(self.a - other.a, self.b - other.b)

    def __neg__(self):
        return Oct(-self.a, -self.b)

    def __mul__(self, other):
        if isinstance(other, Zphi):
            return Oct(self.a * other, self.b * other)
        # Cayley-Dickson: (a, b)(c, d) = (ac - conj(d) b, d a + b conj(c))
        a, b = self.a, self.b
        c, d = other.a, other.b
        return Oct(a * c - d.conj() * b, d * a + b * c.conj())

    __rmul__ = __mul__

    def conj(self) -> "Oct":
        return Oct(self.a.conj(), -self.b)

    def real_part(self) -> Zphi:
        # Re((a, b)) = Re(a)
        return self.a.real_part()

    def norm(self) -> Zphi:
        # N((a, b)) = a a^* + b b^* = N_H(a) + N_H(b)
        return self.a.norm() + self.b.norm()

    def is_zero(self) -> bool:
        return self.a.is_zero() and self.b.is_zero()


# Build the G_0 basis (1, i, h, g, ell, i*ell, h*ell, g*ell)
def _g0_basis():
    bs = []
    for q in ICOSIAN_BASIS:
        bs.append(Oct(q, Quat.of(0, 0, 0, 0)))  # (q, 0)
    for q in ICOSIAN_BASIS:
        bs.append(Oct(Quat.of(0, 0, 0, 0), q))  # (0, q) = q * ell
    return tuple(bs)


G0_BASIS = _g0_basis()
G0_NAMES = ICOSIAN_NAMES + tuple(n + "*ell" for n in ICOSIAN_NAMES)


def bilinear_B(x: Oct, y: Oct) -> Zphi:
    """B(x, y) = Re(x * conj(y))."""
    return (x * y.conj()).real_part()


def g0_gram() -> List[List[Zphi]]:
    """Compute the 8x8 Gram matrix of B on the G_0 basis."""
    return [[bilinear_B(x, y) for y in G0_BASIS] for x in G0_BASIS]


def is_half_integer_not_integer(z: Zphi) -> bool:
    """True iff z is in (1/2) O_K but not in O_K."""
    return (
        (2 * z.a).denominator == 1
        and (2 * z.b).denominator == 1
        and not z.is_integral()
    )


def g0_gram_half_integer_count() -> int:
    G = g0_gram()
    return sum(
        1
        for row in G
        for z in row
        if is_half_integer_not_integer(z)
    )


# Hermitian J_3(G_0).  Represent X by (alpha_1, alpha_2, alpha_3, x1, x2, x3)
# with alpha_i in O_K (Zphi) and x_i in G_0 (Oct).
@dataclass(frozen=True)
class J3Elt:
    alpha: Tuple[Zphi, Zphi, Zphi]
    x: Tuple[Oct, Oct, Oct]

    def __add__(self, other):
        return J3Elt(
            tuple(self.alpha[i] + other.alpha[i] for i in range(3)),
            tuple(self.x[i] + other.x[i] for i in range(3)),
        )

    def __sub__(self, other):
        return J3Elt(
            tuple(self.alpha[i] - other.alpha[i] for i in range(3)),
            tuple(self.x[i] - other.x[i] for i in range(3)),
        )

    def __mul__(self, other):
        # Scalar by Zphi
        if isinstance(other, Zphi):
            return J3Elt(
                tuple(self.alpha[i] * other for i in range(3)),
                tuple(self.x[i] * other for i in range(3)),
            )
        return NotImplemented

    __rmul__ = __mul__

    def to_matrix(self) -> List[List[Oct]]:
        """Build the 3x3 octonion matrix representation."""
        a0, a1, a2 = self.alpha
        x1, x2, x3 = self.x
        zero_oct = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
        # Diagonal: alpha_i are scalars; embed as Oct(scalar, 0)
        def embed_scalar(z):
            return Oct(Quat(z, ZERO, ZERO, ZERO), Quat.of(0, 0, 0, 0))
        M = [
            [embed_scalar(a0), x3, x2.conj()],
            [x3.conj(), embed_scalar(a1), x1],
            [x2, x1.conj(), embed_scalar(a2)],
        ]
        return M


def _e_ii(i: int) -> J3Elt:
    alpha = [ZERO, ZERO, ZERO]
    alpha[i] = ONE
    zero_oct = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
    return J3Elt(tuple(alpha), (zero_oct, zero_oct, zero_oct))


def _f_ij(i: int, j: int, g: Oct) -> J3Elt:
    """F_ij(g) = E_ij(g) + E_ji(g^*).  Encode by x_k slot.

    Slot convention from the manuscript:
      x_1 is the (2,3) entry, x_2 is the (1,3), x_3 is the (1,2).
    """
    zero = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
    alpha = (ZERO, ZERO, ZERO)
    # Map (i, j) -> slot index k in {1,2,3}, where k is the unused index.
    pair = frozenset({i, j})
    if pair == frozenset({1, 2}):
        k = 2  # x_3 entry
        return J3Elt(alpha, (zero, zero, g))
    if pair == frozenset({0, 2}):
        k = 1  # x_2 entry
        return J3Elt(alpha, (zero, g, zero))
    if pair == frozenset({0, 1}):
        k = 0  # x_3 entry... actually slot mapping:
        # In manuscript: (1,2) -> x_3, (1,3) -> x_2, (2,3) -> x_1
        # With 0-indexed i,j: (0,1) -> x_3 = index 2, (0,2) -> x_2 = index 1,
        # (1,2) -> x_1 = index 0
        return J3Elt(alpha, (zero, zero, g))
    raise ValueError(f"bad pair ({i},{j})")


def build_j3_basis() -> List[Tuple[str, J3Elt]]:
    """Return the 27-element ordered basis of J_3(G_0) with labels."""
    basis: List[Tuple[str, J3Elt]] = []
    # Three diagonal idempotents
    for i in range(3):
        basis.append((f"E_{i+1}{i+1}", _e_ii(i)))
    # 24 off-diagonal F_ij(e_a)
    zero_oct = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
    # Slot 12 -> x_3 (manuscript x_3)
    # Slot 13 -> x_2
    # Slot 23 -> x_1
    slot_pairs = [((0, 1), 2), ((0, 2), 1), ((1, 2), 0)]
    for (i, j), idx in slot_pairs:
        for a, name in enumerate(G0_NAMES):
            x = [zero_oct, zero_oct, zero_oct]
            x[idx] = G0_BASIS[a]
            basis.append((f"F_{{{i+1}{j+1}}}({name})", J3Elt(
                (ZERO, ZERO, ZERO), tuple(x))))
    return basis


def matmul_octonion(A: List[List[Oct]], B: List[List[Oct]]) -> List[List[Oct]]:
    """Multiply two 3x3 octonion matrices entry-wise, with standard
    left-to-right parenthesisation of matrix multiplication (no
    associativity assumed inside octonion entries)."""
    n = len(A)
    zero_oct = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
    C = [[zero_oct for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = zero_oct
            for k in range(n):
                s = s + A[i][k] * B[k][j]
            C[i][j] = s
    return C


def is_hermitian_3x3(M: List[List[Oct]]) -> bool:
    """Check M_ij = M_ji^* on a 3x3 octonion matrix and diagonal real."""
    for i in range(3):
        for j in range(3):
            if i == j:
                # Diagonal entry should be a scalar (b component zero).
                if not M[i][i].b.is_zero():
                    return False
            else:
                # M_ij = conj(M_ji)
                diff = M[i][j] - M[j][i].conj()
                if not diff.is_zero():
                    return False
    return True


def doubled_product(X: J3Elt, Y: J3Elt) -> J3Elt:
    """Doubled product {X, Y} = XY + YX, computed by direct matrix
    multiplication followed by hermitian projection back into J_3.

    The result is automatically a hermitian element of J_3(G_0) on
    fully algebraic grounds (the Albert symmetrisation preserves
    hermiticity).  We project back to the (alpha, x) coordinates.
    """
    MX = X.to_matrix()
    MY = Y.to_matrix()
    XY = matmul_octonion(MX, MY)
    YX = matmul_octonion(MY, MX)
    # Sum
    n = 3
    zero_oct = Oct(Quat.of(0, 0, 0, 0), Quat.of(0, 0, 0, 0))
    S = [[XY[i][j] + YX[i][j] for j in range(n)] for i in range(n)]
    return _matrix_to_j3(S)


def _matrix_to_j3(M: List[List[Oct]]) -> J3Elt:
    """Convert a hermitian 3x3 octonion matrix to its (alpha, x) form.

    Slot convention: x_1 = M[1][2], x_2 = M[0][2], x_3 = M[0][1].
    Diagonal entries must have b=0 (scalar embedding).
    """
    a0 = M[0][0].a.x0
    a1 = M[1][1].a.x0
    a2 = M[2][2].a.x0
    x1 = M[1][2]
    x2 = M[0][2]
    x3 = M[0][1]
    return J3Elt((a0, a1, a2), (x1, x2, x3))


def _icosian_coefficients(q: Quat) -> Tuple[Zphi, Zphi, Zphi, Zphi]:
    """Express a quaternion in the icosian basis {1, i, h, g}.

    With h = (1 + i + j + k)/2 and g = (-1 + (phi-1) i - phi j)/2, the
    change of basis from quaternion (1,i,j,k) coefficients (x0,x1,x2,x3)
    to icosian coefficients (c1, ci, ch, cg) is:
       ch = 2 x3
       cg = 2 (phi - 1) (x3 - x2)
       ci = x1 - x3 - (2 - phi)(x3 - x2)
       c1 = x0 - x3 + (phi - 1)(x3 - x2)
    """
    x0, x1, x2, x3 = q.x0, q.x1, q.x2, q.x3
    two = Zphi.of(2, 0)
    phi_minus_one = Zphi.of(-1, 1)        # phi - 1 = phi^{-1}
    two_minus_phi = Zphi.of(2, -1)        # 2 - phi = phi^{-2}
    diff = x3 - x2
    ch = two * x3
    cg = two * phi_minus_one * diff
    ci = x1 - x3 - two_minus_phi * diff
    c1 = x0 - x3 + phi_minus_one * diff
    return (c1, ci, ch, cg)


def _quaternion_in_icosian(q: Quat) -> bool:
    """True iff q is a Z[phi]-linear combination of the icosian basis."""
    return all(c.is_integral() for c in _icosian_coefficients(q))


def _octonion_in_g0(o: Oct) -> bool:
    """True iff (a, b) lies in G_0 = I + I*ell, i.e. both a and b are
    in the icosian ring I."""
    return _quaternion_in_icosian(o.a) and _quaternion_in_icosian(o.b)


def j3_is_integral(X: J3Elt) -> bool:
    """True iff X lies in J_3(G_0): alpha_i in Z[phi] and x_i in G_0."""
    for ai in X.alpha:
        if not ai.is_integral():
            return False
    for xi in X.x:
        if not _octonion_in_g0(xi):
            return False
    return True


def ordinary_jordan_product(X: J3Elt, Y: J3Elt) -> J3Elt:
    """X o Y = (XY + YX) / 2."""
    DP = doubled_product(X, Y)
    return DP * HALF


def jordan_failure_count() -> Tuple[int, int, int, int]:
    """Iterate over the 27x27 = 729 ordered basis pairs and count how
    many fail to be O_K-integral under the ordinary Albert product.

    Returns (DD, DO, OO_same, OO_diff) failure counts.
    """
    basis = build_j3_basis()
    DD = DO = OO_same = OO_diff = 0
    for i, (na, Xa) in enumerate(basis):
        for j, (nb, Xb) in enumerate(basis):
            prod = ordinary_jordan_product(Xa, Xb)
            if not j3_is_integral(prod):
                # Classify by Peirce block
                if i < 3 and j < 3:
                    DD += 1
                elif i < 3 or j < 3:
                    DO += 1
                else:
                    # Both off-diagonal; same-slot or different-slot?
                    # Slots: idx 3..10 = F_12, 11..18 = F_13, 19..26 = F_23
                    def slot(idx):
                        return (idx - 3) // 8
                    if slot(i) == slot(j):
                        OO_same += 1
                    else:
                        OO_diff += 1
    return DD, DO, OO_same, OO_diff


def doubled_failure_count() -> int:
    """Verify that the doubled product closes integrally: should be 0."""
    basis = build_j3_basis()
    failures = 0
    for _, Xa in basis:
        for _, Xb in basis:
            prod = doubled_product(Xa, Xb)
            if not j3_is_integral(prod):
                failures += 1
    return failures


if __name__ == "__main__":
    # Smoke test: G_0 Gram half-integer count
    half_int = g0_gram_half_integer_count()
    print(f"G_0 Gram: half-integer-not-integer entry count = {half_int}")
    print(f"Expected 20 per Lemma 4.4.  {'PASS' if half_int == 20 else 'FAIL'}")
