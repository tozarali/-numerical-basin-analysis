"""Exact linear amplification analysis of the alpha/Newmark family.

Linear test problem: x'' + 2 xi W x' + W^2 x = 0, step h, Omega = W h.
State (x, v, a) -> (x1, v1, a1) through the same equations as
integrators.step, so results apply directly to `parameters(...)` tuples.
"""
import numpy as np


def amplification_matrix(p, omega_h, xi=0.0):
    """3x3 map for unit frequency (W=1) and step omega_h."""
    if not np.isfinite(omega_h) or omega_h < 0 or not np.isfinite(xi) or xi < 0:
        raise ValueError("omega_h and xi must be finite and nonnegative")
    am, af, gamma, beta = p
    h, c, k = float(omega_h), 2.0 * xi, 1.0
    # Unknowns z1=(x1,v1,a1); A z1 = B z0
    A = np.array([[1.0, 0.0, -beta*h*h],
                  [0.0, 1.0, -gamma*h],
                  [(1-af)*k, (1-af)*c, 1-am]])
    B = np.array([[1.0, h, (0.5-beta)*h*h],
                  [0.0, 1.0, (1-gamma)*h],
                  [-af*k, -af*c, -am]])
    return np.linalg.solve(A, B)


def principal_pair(p, omega_h, xi=0.0):
    """Return the principal (oscillatory) eigenvalue with positive argument,
    or None if the principal roots are real (overdamped numerically)."""
    matrix = amplification_matrix(p, omega_h, xi)
    if xi >= 1 or omega_h == 0:
        return None
    lam = np.linalg.eigvals(matrix)
    complex_roots = lam[np.abs(lam.imag) > 1e-12]
    if complex_roots.size == 0:
        return None
    return complex_roots[np.argmax(complex_roots.imag)]


def spectral_radius(p, omega_h, xi=0.0):
    return float(np.max(np.abs(np.linalg.eigvals(amplification_matrix(p, omega_h, xi)))))


def dissipation_and_period(p, omega_h, xi=0.0):
    """Algorithmic damping ratio xi_bar (total, includes physical xi if given)
    and relative period error dT/T = Omega/Omega_bar - 1.
    lambda = exp(Omega_bar h (-xi_bar + i sqrt(1-xi_bar^2)))."""
    lam = principal_pair(p, omega_h, xi)
    if lam is None:
        return {"xi_bar": np.nan, "period_error": np.nan, "modulus": np.nan}
    r, theta = abs(lam), np.angle(lam)
    s = -np.log(r) / theta                    # = xi_bar / sqrt(1-xi_bar^2)
    xi_bar = s / np.sqrt(1 + s*s)
    omega_bar_h = theta / np.sqrt(1 - xi_bar**2)
    return {"xi_bar": float(xi_bar), "period_error": float(omega_h / omega_bar_h - 1),
            "modulus": float(r)}


def per_cycle_energy_ratio(p, omega_h, xi=0.0):
    """Energy retained after one physical period 2*pi (W=1): |lambda|^(2N), N=2pi/h."""
    lam = principal_pair(p, omega_h, xi)
    if lam is None:
        return np.nan
    return float(abs(lam) ** (2 * 2*np.pi / omega_h))


def curves(p, omega_h_values, xi=0.0):
    rows = [dict(omega_h=float(w), **dissipation_and_period(p, w, xi),
                 rho=spectral_radius(p, w, xi)) for w in omega_h_values]
    return rows
