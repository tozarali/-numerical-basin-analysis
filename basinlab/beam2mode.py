"""Two-mode Galerkin model of an axially loaded, simply supported beam with
mid-plane stretching (nondimensional). Modes sin(j*pi*s), j in `modes`.
V(q) = sum_j 1/2 j^2 (j^2 - P) q_j^2 + alpha/4 (sum_j j^2 q_j^2)^2
P in (1, 4): mode 1 buckled (bistable), higher modes stable.
Coupling strength between modes is set by alpha and P (check with the
coupling diagnostics already used for the truss)."""
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Beam2Mode:
    P: float = 2.0
    alpha: float = 1.0
    modes: tuple = (1, 2)
    damping: tuple = (0.15, 0.15)   # viscous per mode coordinate

    def __post_init__(self):
        if not 1 < self.P < self.modes[1]**2:
            raise ValueError("need 1 < P < j2^2 so only mode 1 buckles")

    @property
    def j2(self):
        return np.asarray(self.modes, float)**2

    def potential(self, q):
        j2 = self.j2
        lin = 0.5*np.sum(j2*(j2 - self.P)*q*q, axis=-1)
        s = np.sum(j2*q*q, axis=-1)
        return lin + 0.25*self.alpha*s*s

    def force_tangent(self, q):
        """Internal force dV/dq and Hessian, batched over leading axes."""
        j2 = self.j2
        s = np.sum(j2*q*q, axis=-1)[..., None]
        f = j2*(j2 - self.P)*q + self.alpha*s*j2*q
        K = (np.eye(len(j2))*(j2*(j2 - self.P)) + self.alpha*s[..., None]*np.diag(j2)
             + 2*self.alpha*np.einsum("...i,...j->...ij", j2*q, j2*q))
        return f, K

    def equilibria(self):
        """Buckled states on mode 1: q1 = +/- sqrt((P-1)/alpha), q2 = 0."""
        q1 = np.sqrt((self.P - 1)/self.alpha)
        return np.array([[q1, 0.0], [-q1, 0.0]])

    @property
    def barrier(self):
        """Saddle at q=0 (flat state), V=0; well depth = -(P-1)^2/(4 alpha)."""
        return 0.0

    def energy(self, q, v):
        return self.potential(q) + 0.5*np.sum(v*v, axis=-1)

    def classify(self, q, v):
        e = self.energy(q, v)
        out = np.zeros(e.shape, dtype=np.int8)
        trapped = np.isfinite(e) & (e < self.barrier)
        out[trapped & (q[..., 0] > 0)] = 1
        out[trapped & (q[..., 0] < 0)] = -1
        return out
