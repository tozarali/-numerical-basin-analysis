"""Tilted double well as an MDOF-style model (d = 1) so that the generic
integrators (mdof_genalpha, reference_mdof, dg_mdof, bathe) apply unchanged.

V(x) = (x^2-1)^2/4 - eps*x,  force f = V' = x^3 - x - eps.
eps = 0 is the symmetric well used everywhere else (barrier 1/4).
For eps != 0 the saddle is the middle real root x_s of f; trajectories with
E = v^2/2 + V(x) < V(x_s) are trapped and labelled by sign(x - x_s)."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Well:
    eps: float = 0.0

    def __post_init__(self):
        if abs(self.eps) >= 2/(3*np.sqrt(3)):     # beyond this the wells merge
            raise ValueError("|eps| must be < 2/(3 sqrt 3) for two wells")

    def _V(self, x):
        return 0.25*(x*x - 1)**2 - self.eps*x

    def potential(self, q):
        return self._V(np.asarray(q, float)[..., 0])

    def force(self, q):
        q = np.asarray(q, float)
        return q**3 - q - self.eps

    def force_tangent(self, q):
        q = np.asarray(q, float)
        return q**3 - q - self.eps, (3*q*q - 1)[..., None]

    @property
    def equilibria(self):
        r = np.sort(np.roots([1.0, 0.0, -1.0, -self.eps]).real)
        return r                                  # [left well, saddle, right well]

    @property
    def saddle(self):
        return float(self.equilibria[1])

    @property
    def barrier(self):
        return float(self._V(self.saddle))

    def energy(self, q, v):
        return self.potential(q) + 0.5*np.sum(np.asarray(v, float)**2, axis=-1)

    def classify(self, q, v):
        q, v = np.asarray(q, float), np.asarray(v, float)
        e = self.energy(q, v)
        out = np.zeros(e.shape, dtype=np.int8)
        ok = np.isfinite(e) & (e < self.barrier)
        x = q[..., 0]
        out[ok & (x > self.saddle)] = 1
        out[ok & (x < self.saddle)] = -1
        return out
