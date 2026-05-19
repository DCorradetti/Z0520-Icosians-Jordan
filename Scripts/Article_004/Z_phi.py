"""
Z_phi.py — Exact arithmetic in Z[phi] and Q[phi].

Elements are pairs (a, b) representing a + b*phi, with phi^2 = phi + 1.
Coefficients a, b are exact rationals via fractions.Fraction.

This module supports addition, subtraction, multiplication, scalar
multiplication, Galois conjugation phi -> 1 - phi, integrality check,
numerical evaluation at the two real embeddings sigma_+ and sigma_-,
and hashing for use as dictionary keys.
"""

from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import sqrt


def _frac(x) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    return Fraction(x)


@dataclass(frozen=True)
class Zphi:
    """Element a + b*phi of Q(phi); integral in Z[phi] iff a, b in Z."""
    a: Fraction
    b: Fraction

    @classmethod
    def of(cls, a, b) -> "Zphi":
        return cls(_frac(a), _frac(b))

    def __add__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Zphi.of(other, 0)
        return Zphi(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Zphi.of(other, 0)
        return Zphi(self.a - other.a, self.b - other.b)

    def __rsub__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Zphi.of(other, 0)
        return other - self

    def __neg__(self):
        return Zphi(-self.a, -self.b)

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Zphi.of(other, 0)
        # (a1 + b1 phi)(a2 + b2 phi)
        #   = a1 a2 + (a1 b2 + b1 a2) phi + b1 b2 phi^2
        #   = a1 a2 + b1 b2 + (a1 b2 + b1 a2 + b1 b2) phi
        return Zphi(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Zphi.of(other, 0)
        # 1 / (a + b phi) when phi^2 = phi + 1:
        #   norm = a^2 + a b - b^2 = (a + b phi)(a + b (1 - phi)) reduced
        # Use Galois: conj = a + b - b phi (since phi -> 1 - phi)
        c = other.galois()
        denom = (other * c).a  # b-coefficient is zero
        return self * Zphi(c.a / denom, c.b / denom)

    def galois(self) -> "Zphi":
        """Galois conjugation phi -> 1 - phi."""
        return Zphi(self.a + self.b, -self.b)

    def is_integral(self) -> bool:
        return self.a.denominator == 1 and self.b.denominator == 1

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def sigma_plus(self) -> float:
        return float(self.a) + float(self.b) * (1 + sqrt(5)) / 2

    def sigma_minus(self) -> float:
        return float(self.a) + float(self.b) * (1 - sqrt(5)) / 2

    def __repr__(self):
        if self.b == 0:
            return f"{self.a}"
        if self.a == 0:
            return f"{self.b}*phi"
        return f"({self.a} + {self.b}*phi)"


# Constants
ZERO = Zphi.of(0, 0)
ONE = Zphi.of(1, 0)
PHI = Zphi.of(0, 1)
HALF = Zphi.of(Fraction(1, 2), 0)
PHI_INV = Zphi.of(-1, 1)        # phi^(-1) = phi - 1
TWO = Zphi.of(2, 0)
TWO_MINUS_PHI = Zphi.of(2, -1)   # equals phi^(-2)
PHI_SQ = Zphi.of(1, 1)           # phi^2 = phi + 1


def smoke_test():
    """Minimal smoke test for the arithmetic."""
    # phi^2 = phi + 1
    assert PHI * PHI == PHI_SQ, f"phi^2 = {PHI*PHI}, expected {PHI_SQ}"
    # phi * (phi - 1) = 1
    assert PHI * PHI_INV == ONE
    # 2 - phi = phi^(-2)
    assert PHI_INV * PHI_INV == TWO_MINUS_PHI
    # Galois of phi is 1 - phi
    assert PHI.galois() == Zphi.of(1, -1)
    # Norm = a*(a+b) - b^2*0 + ... check that x * x.galois has b=0
    x = Zphi.of(3, 7)
    norm = x * x.galois()
    assert norm.b == 0, f"norm got b coefficient: {norm}"
    print("Z_phi smoke test PASSED")


if __name__ == "__main__":
    smoke_test()
