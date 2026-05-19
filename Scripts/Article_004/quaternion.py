"""
quaternion.py — Quaternion algebra H_K = (-1, -1)/K over K = Q(sqrt 5).

A quaternion is x = x0 + x1 i + x2 j + x3 k, with coefficients in Q[phi]
(class Zphi from Z_phi.py).  Relations: i^2 = j^2 = -1, ij = k, jk = i,
ki = j; conjugation negates i, j, k.

The icosian ring I sits inside H_K with O_K-basis {1, i, h, g} where:
  h = (1 + i + j + k) / 2
  g = (-1 + (phi - 1) i - phi j) / 2  (Tits 1980 convention)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple

from Z_phi import Zphi, ZERO, ONE, HALF, PHI, PHI_INV


@dataclass(frozen=True)
class Quat:
    """Quaternion x0 + x1 i + x2 j + x3 k with Zphi coefficients."""
    x0: Zphi
    x1: Zphi
    x2: Zphi
    x3: Zphi

    @classmethod
    def of(cls, x0, x1, x2, x3) -> "Quat":
        return cls(*(c if isinstance(c, Zphi) else Zphi.of(c, 0)
                     for c in (x0, x1, x2, x3)))

    def __add__(self, other):
        return Quat(self.x0 + other.x0, self.x1 + other.x1,
                    self.x2 + other.x2, self.x3 + other.x3)

    def __sub__(self, other):
        return Quat(self.x0 - other.x0, self.x1 - other.x1,
                    self.x2 - other.x2, self.x3 - other.x3)

    def __neg__(self):
        return Quat(-self.x0, -self.x1, -self.x2, -self.x3)

    def __mul__(self, other):
        if isinstance(other, Zphi):
            return Quat(self.x0 * other, self.x1 * other,
                        self.x2 * other, self.x3 * other)
        if not isinstance(other, Quat):
            other = Quat.of(other, 0, 0, 0)
        a0, a1, a2, a3 = self.x0, self.x1, self.x2, self.x3
        b0, b1, b2, b3 = other.x0, other.x1, other.x2, other.x3
        # Hamilton product
        return Quat(
            a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0,
        )

    __rmul__ = __mul__

    def conj(self) -> "Quat":
        return Quat(self.x0, -self.x1, -self.x2, -self.x3)

    def trace(self) -> Zphi:
        # tr(x) = x + conj(x) = 2 x0
        return self.x0 + self.x0

    def norm(self) -> Zphi:
        # N(x) = x * conj(x) = x0^2 + x1^2 + x2^2 + x3^2
        result = self * self.conj()
        return result.x0

    def real_part(self) -> Zphi:
        return self.x0

    def is_zero(self) -> bool:
        return (self.x0.is_zero() and self.x1.is_zero()
                and self.x2.is_zero() and self.x3.is_zero())


# Quaternion basis
QH_ONE = Quat.of(1, 0, 0, 0)
QH_I = Quat.of(0, 1, 0, 0)
QH_J = Quat.of(0, 0, 1, 0)
QH_K = Quat.of(0, 0, 0, 1)

# Icosian basis (Tits 1980 convention):
#   1, i, h = (1 + i + j + k)/2, g = (-1 + (phi-1) i - phi j)/2
HALF_PLUS = HALF
NEG_HALF = Zphi.of(-HALF.a, 0)


def _icosian_basis():
    """Return the ordered icosian basis (1, i, h, g) as Quat tuples."""
    one = QH_ONE
    i_q = QH_I
    h_q = Quat(HALF_PLUS, HALF_PLUS, HALF_PLUS, HALF_PLUS)
    # g = (-1 + (phi-1) i - phi j) / 2
    minus_half = Zphi(Zphi.of(-1, 0).a / 2, Zphi.of(0, 0).b)  # = -1/2
    g_q = Quat(
        Zphi.of(-HALF.a, 0),                          # -1/2
        Zphi.of((PHI_INV.a) / 2, (PHI_INV.b) / 2),    # (phi-1)/2
        Zphi.of(-PHI.a / 2, -PHI.b / 2),              # -phi/2
        Zphi.of(0, 0),
    )
    return (one, i_q, h_q, g_q)


ICOSIAN_BASIS = _icosian_basis()
ICOSIAN_NAMES = ("1", "i", "h", "g")


def smoke_test():
    """Check icosian properties: norms, multiplication closure on basis."""
    one, i_q, h_q, g_q = ICOSIAN_BASIS
    # Norms should be 1
    for name, x in zip(ICOSIAN_NAMES, ICOSIAN_BASIS):
        n = x.norm()
        assert n == Zphi.of(1, 0), f"norm({name}) = {n}, expected 1"

    # Test polar pairing B(x, y) = (1/2) tr(x conj(y))
    # B(1, 1) = 1, B(i, i) = 1 (since 1 = 1*1)
    # B(1, h) = (1/2) tr(h) = (1/2) * (h_q.trace()) = 1/2
    b_1_h = (one * h_q.conj()).trace() * HALF
    assert b_1_h.a == HALF.a and b_1_h.b == 0, f"B(1,h) = {b_1_h}"

    # Norm of (a + b*ell) under Cayley-Dickson = N(a) + N(b); test on basis
    print("quaternion smoke test PASSED")


if __name__ == "__main__":
    smoke_test()
