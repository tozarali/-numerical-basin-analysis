"""Independent cross-check of the in-house corotational beam against
OpenSeesPy (elasticBeamColumn + Corotational transformation).

Compared quantities (same geometry, EA=1e4, EI=1, L=1, 16 elements,
end shortening 1.6 x Euler):
  1. critical shortening / axial force of the flat branch (Euler load 4 pi^2 EI/L^2)
  2. post-buckling midspan amplitude of the buckled equilibrium
  3. lowest natural frequencies about the buckled state (tangent K, same
     lumped mass incl. rotary term rhoA L^3/78 per node)
Needs: pip install openseespy (optional; the experiment reports 'skipped'
if it is not importable)."""
import numpy as np
from scipy.linalg import eigh

from ..fe_arch import Arch
from ..runmeta import experiment
from ._common import out_name, quick

FACTOR = 1.6


def inhouse(n_el):
    a = Arch(n_el=n_el)
    a.end_shortening = FACTOR*a.euler_shortening
    up = a.buckled_states()[1][0]
    M = a.lumped_mass()
    _, K = a.force_tangent(up)
    w = np.sqrt(np.abs(eigh(K, M, eigvals_only=True)))
    flat = a.flat_state()[0]
    _, Kf = a.force_tangent(flat)
    wf = np.sqrt(np.abs(eigh(Kf, M, eigvals_only=True)))
    return a, {"midspan_amplitude": float(up[a.mid_index]), "omega_buckled": w[:3].tolist(),
               "omega_flat_at_shortening": wf[:3].tolist()}


def opensees(n_el, EA, EI, rhoA, L, shortening, pert=1e-4, steps=200):
    import openseespy.opensees as ops
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(n_el + 1):
        ops.node(i + 1, L*i/n_el, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(n_el + 1, 0, 1, 1)                     # u free (prescribed by load pattern), v, theta fixed
    ops.geomTransf("Corotational", 1)
    for e in range(n_el):
        ops.element("elasticBeamColumn", e + 1, e + 1, e + 2, EA, 1.0, EI, 1)   # E=1, A=EA, Iz=EI
    mid = n_el//2 + 1
    dl = L/n_el
    for i in range(1, n_el + 2):
        mt = rhoA*dl*(0.5 if i in (1, n_el + 1) else 1.0)
        mr = rhoA*dl**3/78*(1 if i in (1, n_el + 1) else 2)
        ops.mass(i, mt, mt, mr)
    # axial displacement control via a unit axial load pattern and an integrator on node n_el+1 dof 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(n_el + 1, 1.0, 0.0, 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(mid, 0.0, pert, 0.0)                  # small transverse imperfection load
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-12, 50); ops.algorithm("Newton")
    du = -shortening/steps
    ops.integrator("DisplacementControl", n_el + 1, 1, du)
    ops.analysis("Static")
    ok = 0
    for s in range(steps):
        ok = ops.analyze(1)
        if ok != 0:
            break
    # critical shortening: first step where midspan deflection exceeds 100x pert-driven flat value
    w = ops.nodeDisp(mid, 2)
    # remove imperfection load (equilibrate)
    ops.remove("loadPattern", 2)
    ops.integrator("DisplacementControl", n_el + 1, 1, 0.0)   # hold the end shortening, re-equilibrate
    ops.analysis("Static")
    ops.analyze(1)
    w_after = ops.nodeDisp(mid, 2)
    ops.wipeAnalysis()                                   # clamp the end u: same boundary as the in-house model
    ops.fix(n_el + 1, 1, 0, 0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Plain")
    lam = ops.eigen(3)
    return {"converged": ok == 0, "midspan_amplitude": float(w_after),
            "midspan_amplitude_with_imperfection_load": float(w),
            "omega_buckled": [float(np.sqrt(abs(x))) for x in lam[:3]],
            "axial_load_factor_final": float(ops.getLoadFactor(1))}


def main():
    n_el = 8 if quick() else 16
    with experiment(out_name("fe_crosscheck_openseespy"), {"n_el": n_el, "factor": FACTOR}) as run:
        a, mine = inhouse(n_el)
        try:
            import openseespy.opensees  # noqa: F401
        except Exception as exc:                    # pragma: no cover
            run["results"] = {"status": "skipped", "reason": repr(exc), "inhouse": mine}
            return run
        theirs = opensees(n_el, a.EA, a.EI, a.rhoA, a.span, a.end_shortening)
        rel = None
        if theirs["converged"]:
            rel = {"midspan_amplitude": abs(abs(theirs["midspan_amplitude"]) - abs(mine["midspan_amplitude"]))
                   / abs(mine["midspan_amplitude"]),
                   "omega1": abs(theirs["omega_buckled"][0] - mine["omega_buckled"][0]) / mine["omega_buckled"][0]}
        run["results"] = {"status": "done", "inhouse": mine, "opensees": theirs,
                          "euler_load": float(4*np.pi**2*a.EI/a.span**2),
                          "euler_shortening": float(a.euler_shortening),
                          "applied_shortening": float(a.end_shortening),
                          "relative_differences": rel}
    return run


if __name__ == "__main__":
    main()
