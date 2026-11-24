"""Generalized-alpha family, with explicit Newton convergence checks."""
import numpy as np

from .model import force, tangent


def parameters(family, value):
    if family == "generalized-alpha":
        if not 0 <= value <= 1:
            raise ValueError("rho must be in [0, 1]")
        am, af = (2 * value - 1) / (value + 1), value / (value + 1)
    elif family == "HHT":
        if not 0 <= value <= 1 / 3:
            raise ValueError("HHT parameter must be in [0, 1/3]")
        am, af = 0.0, value
    elif family == "WBZ":
        if not -1 <= value <= 0:
            raise ValueError("WBZ parameter must be in [-1, 0]")
        am, af = value, 0.0
    elif family == "Newmark":
        if value < 0:
            raise ValueError("Newmark parameter must be nonnegative")
        return 0.0, 0.0, 0.5 + value, 0.25 * (1 + value)**2
    else:
        raise ValueError(f"Unknown family: {family}")
    return am, af, 0.5 - am + af, 0.25 * (1 - am + af)**2


def spectral_radius_infinity(p):
    """Measure the scheme's high-frequency spectral radius by its linear map."""
    h = 1e3
    am, af, gamma, beta = p
    c1 = 1 / (beta * h * h)
    matrix = np.zeros((3, 3))
    for column in range(3):
        x, v, a = np.eye(3)[:, column]
        x1 = x + h * v + 0.5 * h * h * a
        for _ in range(5):
            a1 = c1 * (x1 - x - h * v) - (0.5 / beta - 1) * a
            v1 = v + h * ((1 - gamma) * a + gamma * a1)
            residual = ((1 - am) * a1 + am * a
                        + (1 - af) * x1 + af * x)
            if abs(residual) <= 1e-12:
                break
            tangent_value = (1 - am) * c1 + (1 - af)
            x1 -= residual / tangent_value
        else:
            raise RuntimeError("Could not construct high-frequency amplification map")
        a1 = c1 * (x1 - x - h * v) - (0.5 / beta - 1) * a
        v1 = v + h * ((1 - gamma) * a + gamma * a1)
        matrix[:, column] = x1, v1, a1
    return float(np.max(np.abs(np.linalg.eigvals(matrix))))


def parameters_for_spectral_radius(family, target):
    if not 0 <= target <= 1:
        raise ValueError("target spectral radius must be in [0, 1]")
    if family == "generalized-alpha":
        return parameters(family, target), target
    ranges = {"HHT": (0.0, 1 / 3), "WBZ": (-1.0, 0.0),
              "Newmark": (0.0, 3.0)}
    if family not in ranges:
        raise ValueError(f"Unknown family: {family}")
    low, high = ranges[family]
    rho_low = spectral_radius_infinity(parameters(family, low))
    rho_high = spectral_radius_infinity(parameters(family, high))
    if not min(rho_low, rho_high) - 1e-9 <= target <= max(rho_low, rho_high) + 1e-9:
        raise ValueError(f"{family} cannot attain spectral radius {target}")
    for _ in range(60):
        middle = (low + high) / 2
        rho_middle = spectral_radius_infinity(parameters(family, middle))
        if (rho_middle - target) * (rho_low - target) > 0:
            low, rho_low = middle, rho_middle
        else:
            high, rho_high = middle, rho_middle
    value = (low + high) / 2
    return parameters(family, value), spectral_radius_infinity(parameters(family, value))


def step(x, v, a, h, damping, p, tol=1e-10, maxiter=30):
    am, af, gamma, beta = p
    c1 = 1 / (beta * h * h)
    x1 = x + h * v + 0.5 * h * h * a
    f0 = force(x)
    for iteration in range(maxiter + 1):
        a1 = c1 * (x1 - x - h * v) - (0.5 / beta - 1) * a
        v1 = v + h * ((1 - gamma) * a + gamma * a1)
        residual = ((1 - am) * a1 + am * a
                    + damping * ((1 - af) * v1 + af * v)
                    + (1 - af) * force(x1) + af * f0)
        if not np.all(np.isfinite(residual)):
            raise FloatingPointError("Non-finite generalized-alpha residual")
        active = np.abs(residual) > tol
        if not np.any(active):
            return x1, v1, a1
        if iteration == maxiter:
            break
        derivative = ((1 - am) * c1 + damping * (1 - af) * gamma * h * c1
                      + (1 - af) * tangent(x1))
        if np.any(active & (np.abs(derivative) < 1e-14)):
            raise ArithmeticError("Singular generalized-alpha Newton tangent")
        x1 = np.where(active, x1 - residual / derivative, x1)
    raise RuntimeError(f"Generalized-alpha Newton iteration failed after {maxiter} steps")


def integrate(x, v, end, h, damping, p):
    if h <= 0 or end <= 0 or damping <= 0:
        raise ValueError("end, h and damping must be positive")
    steps = round(end / h)
    if steps < 1 or not np.isclose(steps * h, end, rtol=0, atol=1e-10):
        raise ValueError("end must be an integer multiple of h")
    x, v = np.asarray(x, dtype=float).copy(), np.asarray(v, dtype=float).copy()
    a = -damping * v - force(x)
    for _ in range(steps):
        x, v, a = step(x, v, a, h, damping, p)
    return x, v
