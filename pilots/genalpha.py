"""Vectorised SDOF generalized-alpha (Chung-Hulbert) with an explicit inner Newton loop.
m a + d v + f(x) = F(t). Inner-loop settings (tol, maxit, tangent) are exposed on purpose."""
import numpy as np

def params(rho):
    am = (2 * rho - 1) / (rho + 1); af = rho / (rho + 1)
    g = 0.5 - am + af; b = 0.25 * (1 - am + af) ** 2
    return am, af, g, b

def step(x0, v0, a0, t0, h, f, k, d=0.0, F=lambda t: 0.0, rho=1.0,
         tol=1e-12, maxit=30, tangent="full", kfix=None):
    am, af, g, b = params(rho)
    c1 = 1.0 / (b * h * h)
    x1 = x0 + h * v0 + 0.5 * h * h * a0          # explicit predictor
    Ft = F(t0 + (1 - af) * h); f0 = f(x0)
    if tangent == "step": kfix = k(x0)
    it = np.zeros(np.shape(x0), int)
    for _ in range(maxit):
        a1 = c1 * (x1 - x0 - h * v0) - (0.5 / b - 1) * a0
        v1 = v0 + h * ((1 - g) * a0 + g * a1)
        R = (1 - am) * a1 + am * a0 + d * ((1 - af) * v1 + af * v0) + (1 - af) * f(x1) + af * f0 - Ft
        act = np.abs(R) > tol
        if not np.any(act): break
        kt = k(x1) if tangent == "full" else kfix
        dR = (1 - am) * c1 + d * (1 - af) * g * h * c1 + (1 - af) * kt
        x1 = np.where(act, x1 - R / dR, x1); it = it + act
    a1 = c1 * (x1 - x0 - h * v0) - (0.5 / b - 1) * a0
    v1 = v0 + h * ((1 - g) * a0 + g * a1)
    return x1, v1, a1, it
