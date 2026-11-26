"""Shared runner for the extended-controls experiments: one call integrates
a batch of ICs with a named method and returns labels."""
import numpy as np

from .. import bathe, dg_mdof
from ..mdof_genalpha import integrate_to_masked as ga_integrate
from ._common import family_parameters


def method_list(rho_step=0.1):
    """(name, kind, rho) tuples: tuned generalized-alpha, rho-Bathe, and
    the two parameter-free controls (trapezoid = gen-alpha rho=1, DG)."""
    rhos = [round(float(r), 4) for r in np.arange(0.0, 1.0 + 1e-9, rho_step)]
    out = [("generalized-alpha", "ga", r) for r in rhos]
    out += [("rho-Bathe", "bathe", r) for r in rhos]
    out += [("discrete-gradient", "dg", None)]
    return out


def integrate(kind, rho, q0, v0, T, h, M, C, model):
    if kind == "ga":
        p, _ = family_parameters("generalized-alpha", rho)
        return ga_integrate(q0, v0, T, h, M, C, model, p)
    if kind == "bathe":
        return bathe.integrate_to_masked(q0, v0, T, h, M, C, model, rho)
    if kind == "dg":
        return dg_mdof.integrate_to_masked(q0, v0, T, h, M, C, model)
    raise ValueError(kind)


def key_of(kind, rho):
    return f"{kind}_{rho}"


def summarize(rows, group_keys):
    """Per group: best generalized-alpha, best Bathe, trapezoid (ga rho=1),
    DG, using the valid-count-weighted mismatch fractions stored in rows."""
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[k] for k in group_keys)].append(r)
    out = []
    for k, rs in g.items():
        def best(kind):
            c = [r for r in rs if r["kind"] == kind and r["mismatch"] is not None]
            return min(c, key=lambda r: r["mismatch"]) if c else None
        def at(kind, rho):
            c = [r for r in rs if r["kind"] == kind and r.get("rho") == rho and r["mismatch"] is not None]
            return c[0]["mismatch"] if c else None
        bg, bb = best("ga"), best("bathe")
        dg = [r for r in rs if r["kind"] == "dg" and r["mismatch"] is not None]
        out.append(dict(zip(group_keys, k),
                        trapezoid=at("ga", 1.0), dg=dg[0]["mismatch"] if dg else None,
                        best_ga=bg["mismatch"] if bg else None, best_ga_rho=bg["rho"] if bg else None,
                        best_bathe=bb["mismatch"] if bb else None, best_bathe_rho=bb["rho"] if bb else None,
                        bathe_rho0=at("bathe", 0.0), ga_rho0=at("ga", 0.0)))
    return out
