"""Geometrically nonlinear, two-bar von Mises truss with a unit apex mass."""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .integrators import parameters


@dataclass(frozen=True)
class Truss:
    half_span: float = 1.0
    rise: float = 0.7
    axial_stiffness: float = 1.0
    damping: float = 0.15

    def __post_init__(self):
        if not (0 < self.rise < np.sqrt(3) * self.half_span
                and self.half_span > 0 and self.axial_stiffness > 0
                and self.damping > 0):
            raise ValueError("Require positive parameters and rise < sqrt(3)*half_span")

    @property
    def rest_length(self):
        return float(np.hypot(self.half_span, self.rise))

    @property
    def barrier(self):
        return self.axial_stiffness * (self.rest_length - self.half_span)**2

    def potential(self, q):
        u, y = q[..., 0], q[..., 1]
        left = np.hypot(u + self.half_span, y)
        right = np.hypot(u - self.half_span, y)
        return 0.5 * self.axial_stiffness * (
            (left - self.rest_length)**2 + (right - self.rest_length)**2)

    def force_tangent(self, q):
        u, y = q[..., 0], q[..., 1]
        force = np.zeros_like(q)
        tangent = np.zeros(q.shape[:-1] + (2, 2))
        for offset in (self.half_span, -self.half_span):
            dx = u + offset
            length = np.hypot(dx, y)
            if np.any(length < 1e-10):
                raise ArithmeticError("A truss bar collapsed to zero length")
            ratio = self.rest_length / length
            force[..., 0] += self.axial_stiffness * (1 - ratio) * dx
            force[..., 1] += self.axial_stiffness * (1 - ratio) * y
            scale = self.axial_stiffness * ratio / length**2
            tangent[..., 0, 0] += self.axial_stiffness * (1-ratio) + scale * dx**2
            tangent[..., 1, 1] += self.axial_stiffness * (1-ratio) + scale * y**2
            tangent[..., 0, 1] += scale * dx * y
        tangent[..., 1, 0] = tangent[..., 0, 1]
        return force, tangent

    def energy(self, q, v):
        return self.potential(q) + 0.5 * np.sum(v*v, axis=-1)

    def classify(self, q, v):
        """Return the well sign only after energy drops below the y=0 barrier."""
        energy = self.energy(q, v)
        resolved = np.isfinite(energy) & (energy < self.barrier)
        return np.where(resolved, np.sign(q[..., 1]), 0).astype(np.int8)


def rk4(model, q, v, end, h, constrained=False):
    if h <= 0 or end <= 0:
        raise ValueError("end and h must be positive")
    steps = round(end / h)
    if steps < 1 or not np.isclose(steps * h, end, rtol=0, atol=1e-10):
        raise ValueError("end must be a positive multiple of h")
    q, v = np.array(q, dtype=float, copy=True), np.array(v, dtype=float, copy=True)

    def rhs(position, velocity):
        force, _ = model.force_tangent(position)
        acceleration = -force - model.damping * velocity
        if constrained:
            acceleration[..., 0] = 0
        return velocity, acceleration

    if constrained and (np.any(q[..., 0] != 0) or np.any(v[..., 0] != 0)):
        raise ValueError("constrained model requires zero horizontal position and velocity")
    for _ in range(steps):
        k1q, k1v = rhs(q, v)
        k2q, k2v = rhs(q+h*k1q/2, v+h*k1v/2)
        k3q, k3v = rhs(q+h*k2q/2, v+h*k2v/2)
        k4q, k4v = rhs(q+h*k3q, v+h*k3v)
        q += h*(k1q+2*k2q+2*k3q+k4q)/6
        v += h*(k1v+2*k2v+2*k3v+k4v)/6
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(v)):
        raise FloatingPointError("Non-finite truss RK4 state")
    return q, v


def generalized_alpha(model, q, v, end, h, rho, constrained=False,
                      tol=1e-10, maxiter=30):
    if h <= 0 or end <= 0:
        raise ValueError("end and h must be positive")
    steps = round(end / h)
    if steps < 1 or not np.isclose(steps * h, end, rtol=0, atol=1e-10):
        raise ValueError("end must be a positive multiple of h")
    q, v = np.array(q, dtype=float, copy=True), np.array(v, dtype=float, copy=True)
    if constrained and (np.any(q[..., 0] != 0) or np.any(v[..., 0] != 0)):
        raise ValueError("constrained model requires zero horizontal position and velocity")
    am, af, gamma, beta = parameters("generalized-alpha", rho)
    c1 = 1 / (beta*h*h)
    f0, _ = model.force_tangent(q)
    a = -f0 - model.damping*v
    if constrained:
        a[..., 0] = 0
    for step_index in range(steps):
        q1 = q + h*v + 0.5*h*h*a
        for iteration in range(maxiter + 1):
            a1 = c1*(q1-q-h*v) - (0.5/beta-1)*a
            v1 = v + h*((1-gamma)*a + gamma*a1)
            f1, kt = model.force_tangent(q1)
            residual = ((1-am)*a1 + am*a
                        + model.damping*((1-af)*v1+af*v)
                        + (1-af)*f1 + af*f0)
            if constrained:
                residual[..., 0] = 0
            if not np.all(np.isfinite(residual)):
                raise FloatingPointError(f"Non-finite truss residual at step {step_index}")
            active = np.max(np.abs(residual), axis=-1) > tol
            if not np.any(active):
                break
            if iteration == maxiter:
                raise RuntimeError(f"Truss Newton failure at step {step_index}: "
                                   f"{int(active.sum())} unconverged points")
            diagonal = (1-am)*c1 + model.damping*(1-af)*gamma*h*c1
            A = (1-af)*kt
            A[..., 0, 0] += diagonal
            A[..., 1, 1] += diagonal
            if constrained:
                correction = np.zeros_like(q)
                correction[..., 1] = residual[..., 1] / A[..., 1, 1]
            else:
                det = A[..., 0, 0]*A[..., 1, 1] - A[..., 0, 1]*A[..., 1, 0]
                if np.any(active & (np.abs(det) < 1e-14)):
                    raise ArithmeticError(f"Singular truss tangent at step {step_index}")
                correction = np.empty_like(q)
                correction[..., 0] = (A[..., 1, 1]*residual[..., 0]
                                      - A[..., 0, 1]*residual[..., 1])/det
                correction[..., 1] = (A[..., 0, 0]*residual[..., 1]
                                      - A[..., 1, 0]*residual[..., 0])/det
            q1 = np.where(active[..., None], q1-correction, q1)
        q, v, a, f0 = q1, v1, a1, f1
    return q, v


def dop853(model, q0, v0, end, rtol=1e-11, atol=1e-13, constrained=False):
    q0, v0 = np.asarray(q0, dtype=float), np.asarray(v0, dtype=float)
    if q0.shape != (2,) or v0.shape != (2,):
        raise ValueError("DOP853 reference takes one 2-component state")
    if constrained and (q0[0] != 0 or v0[0] != 0):
        raise ValueError("constrained model requires zero horizontal initial state")

    def rhs(t, state):
        velocity = state[2:]
        force, _ = model.force_tangent(state[:2])
        acceleration = -force - model.damping*velocity
        if constrained:
            acceleration[0] = 0
        return np.concatenate((velocity, acceleration))

    sol = solve_ivp(rhs, (0, end), np.concatenate((q0, v0)), method="DOP853",
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"DOP853 failed: {sol.message}")
    return sol.y[:2, -1], sol.y[2:, -1]
