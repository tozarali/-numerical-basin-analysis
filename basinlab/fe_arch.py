"""Clamped-clamped beam/arch with 2D corotational Euler-Bernoulli elements
(Crisfield-type), batched over initial conditions. DOFs/node: (u, w, theta).

Why a *buckled* beam: an exploratory check (2026-11-24) showed that a
stress-free shallow clamped sinusoidal arch (rise 0.1, EA/EI=1e4) is not
bistable when unloaded: the midspan reaction never returns to zero along
the displacement-controlled path. Bistability therefore needs residual
compression. Here the right support is moved inwards by `end_shortening`
(prescribed u), which buckles an (optionally imperfect) straight beam into
two mirror-image stable states - the FE analogue of `beam2mode`.

Closed-form checks used in tests/validation (textbook results):
  * clamped-clamped Euler load P_cr = 4 pi^2 EI / L^2, i.e. critical end
    shortening Delta_cr = P_cr L / EA (small strain);
  * first-order post-buckling for w = (w_m/2)(1-cos 2 pi s/L):
    Delta - Delta_cr ~= pi^2 w_m^2 / (4 L).
Independent cross-check against OpenSeesPy: experiments/fe_crosscheck_openseespy.py.
Lumped mass with rotary inertia rhoA*l^3/78 per element end (modelling choice).
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Arch:
    span: float = 1.0
    n_el: int = 16                 # even -> midspan node
    EA: float = 1.0e4
    EI: float = 1.0
    rhoA: float = 1.0
    imperfection: float = 0.0      # initial rise of (1-cos)/2 shape
    end_shortening: float = 0.0    # prescribed inward u at the right end
    rayleigh_mass: float = 0.0     # C = a0 M (+ a1 K0, not used)

    def __post_init__(self):
        if self.n_el % 2:
            raise ValueError("n_el must be even (midspan node)")
        s = np.linspace(0, self.span, self.n_el + 1)
        self.X = np.column_stack((s, self.imperfection*0.5*(1 - np.cos(2*np.pi*s/self.span))))
        self.ndof = 3*(self.n_el + 1)
        self.right_u = self.ndof - 3
        fixed = [0, 1, 2, self.ndof-3, self.ndof-2, self.ndof-1]
        self.free = np.setdiff1d(np.arange(self.ndof), fixed)
        d = self.X[1:] - self.X[:-1]
        self.L0 = np.hypot(d[:, 0], d[:, 1])
        self.beta0 = np.arctan2(d[:, 1], d[:, 0])
        self.mid_w = 3*(self.n_el//2) + 1
        self.mid_index = int(np.flatnonzero(self.free == self.mid_w)[0])

    @property
    def euler_shortening(self):
        return 4*np.pi**2*self.EI/self.span**2*self.span/self.EA

    # --- batched element kinematics: de has shape (B, 6) ---------------------
    def _element(self, e, de):
        X1, X2 = self.X[e], self.X[e+1]
        dx = X2[0] - X1[0] + de[:, 3] - de[:, 0]
        dy = X2[1] - X1[1] + de[:, 4] - de[:, 1]
        Ln = np.hypot(dx, dy)
        beta = np.arctan2(dy, dx)
        rot = np.arctan2(np.sin(beta - self.beta0[e]), np.cos(beta - self.beta0[e]))
        t1, t2 = de[:, 2] - rot, de[:, 5] - rot
        L0 = self.L0[e]
        u = Ln - L0
        N = self.EA*u/L0
        M1 = self.EI/L0*(4*t1 + 2*t2)
        M2 = self.EI/L0*(2*t1 + 4*t2)
        c, s = np.cos(beta), np.sin(beta)
        zero = np.zeros_like(c)
        r = np.stack((-c, -s, zero, c, s, zero), axis=1)
        z = np.stack((s, -c, zero, -s, c, zero), axis=1)
        e3 = np.array([0, 0, 1, 0, 0, 0.0]); e6 = np.array([0, 0, 0, 0, 0, 1.0])
        B = np.stack((r, -z/Ln[:, None] + e3, -z/Ln[:, None] + e6), axis=1)   # (B,3,6)
        return dict(Ln=Ln, u=u, t1=t1, t2=t2, N=N, M1=M1, M2=M2, r=r, z=z, B=B, L0=L0)

    def _assemble(self, d, want_tangent=True, want_energy=False):
        """d: (B, ndof) full displacement vectors."""
        nb = d.shape[0]
        f = np.zeros((nb, self.ndof))
        K = np.zeros((nb, self.ndof, self.ndof)) if want_tangent else None
        E = np.zeros(nb)
        for e in range(self.n_el):
            sl = slice(3*e, 3*e + 6)
            k = self._element(e, d[:, sl])
            sig = np.stack((k["N"], k["M1"], k["M2"]), axis=1)
            f[:, sl] += np.einsum("bij,bi->bj", k["B"], sig)
            if want_tangent:
                L0 = k["L0"]
                D = np.array([[self.EA/L0, 0, 0], [0, 4*self.EI/L0, 2*self.EI/L0],
                              [0, 2*self.EI/L0, 4*self.EI/L0]])
                Ke = np.einsum("bki,kl,blj->bij", k["B"], D, k["B"])
                Ke += (k["N"]/k["Ln"])[:, None, None]*np.einsum("bi,bj->bij", k["z"], k["z"])
                Ke += ((k["M1"] + k["M2"])/k["Ln"]**2)[:, None, None]*(
                    np.einsum("bi,bj->bij", k["r"], k["z"]) + np.einsum("bi,bj->bij", k["z"], k["r"]))
                K[:, sl, sl] += Ke
            if want_energy:
                E += 0.5*(k["N"]*k["u"] + k["M1"]*k["t1"] + k["M2"]*k["t2"])
        return f, K, E

    def internal_force(self, d_full):
        """Full-DOF internal force for one displacement vector (ndof,)."""
        f, _, _ = self._assemble(np.atleast_2d(d_full), want_tangent=False)
        return f[0]

    # --- reduced (free DOF) interface, compatible with mdof_genalpha -------
    def full(self, q):
        q = np.atleast_2d(q)
        d = np.zeros((q.shape[0], self.ndof))
        d[:, self.free] = q
        d[:, self.right_u] = -self.end_shortening
        return d

    def force_tangent(self, q):
        q = np.asarray(q, float)
        shape = q.shape
        f, K, _ = self._assemble(self.full(q.reshape(-1, len(self.free))))
        f = f[:, self.free]
        K = K[:, self.free][:, :, self.free]
        return f.reshape(shape), K.reshape(shape + (len(self.free),))

    def force(self, q):
        q = np.asarray(q, float)
        f, _, _ = self._assemble(self.full(q.reshape(-1, len(self.free))), want_tangent=False)
        return f[:, self.free].reshape(q.shape)

    def potential(self, q):
        q = np.asarray(q, float)
        _, _, E = self._assemble(self.full(q.reshape(-1, len(self.free))),
                                 want_tangent=False, want_energy=True)
        return E.reshape(q.shape[:-1])

    def mass_diag(self):
        m = np.zeros(self.ndof)
        for e, L in enumerate(self.L0):
            for node in (e, e+1):
                m[3*node:3*node+2] += 0.5*self.rhoA*L
                m[3*node+2] += self.rhoA*L**3/78      # rotary term: CHECK/justify
        return m[self.free]

    def lumped_mass(self):
        return np.diag(self.mass_diag())

    def damping_matrix(self):
        return self.rayleigh_mass*self.lumped_mass()

    def energy(self, q, v):
        return self.potential(q) + 0.5*np.sum(self.mass_diag()*v*v, axis=-1)

    # --- statics ------------------------------------------------------------
    def equilibrium(self, q0, tol=1e-9, maxiter=60, load_steps=1):
        """Newton on f(q)=0 (free DOFs). With load_steps>1 the end shortening
        is ramped from 0 to its value (path following for buckling).
        Returns q, converged, sorted eigenvalues of K (smallest first)."""
        q = np.asarray(q0, float).copy()
        target = self.end_shortening
        converged = False
        for j in range(1, load_steps + 1):
            self.end_shortening = target*j/load_steps
            converged = False
            for _ in range(maxiter):
                f, K = self.force_tangent(q)
                if np.linalg.norm(f) < tol*max(1.0, self.EA*self.end_shortening):
                    converged = True
                    break
                q -= np.linalg.solve(K, f)
        self.end_shortening = target
        _, K = self.force_tangent(q)
        return q, converged, np.linalg.eigvalsh(K)

    def shape_vector(self, amplitude=1.0):
        """Free-DOF vector of w = (amplitude/2)(1 - cos 2 pi s/L), theta = w'."""
        s = self.X[:, 0]
        k = 2*np.pi/self.span
        w = 0.5*amplitude*(1 - np.cos(k*s))
        th = 0.5*amplitude*k*np.sin(k*s)
        d = np.zeros(self.ndof)
        d[1::3], d[2::3] = w, th
        return d[self.free]

    def buckled_states(self):
        """Up (+1) and down (-1) states: ramp shortening on the flat branch,
        then start Newton from the first-order post-buckling amplitude."""
        flat, _, _ = self.flat_state()
        excess = max(self.end_shortening - self.euler_shortening, 0.0)
        guess = np.sqrt(4*self.span*excess/np.pi**2) if excess > 0 else 1e-3
        return {sign: self.equilibrium(flat + self.shape_vector(sign*guess))
                for sign in (1, -1)}

    def flat_state(self):
        """Straight compressed equilibrium (symmetric candidate saddle)."""
        return self.equilibrium(np.zeros(len(self.free)), load_steps=5)

    def constrained_midspan(self, w_values, tol=1e-9, maxiter=60):
        """Equilibria with the midspan deflection prescribed (all other free
        DOFs relaxed), by continuation outward from the flat state. Returns
        array (len(w_values), n_free) and the min eigenvalue of the
        constrained tangent (negative -> constrained state unstable)."""
        w_values = np.asarray(w_values, float)
        k = self.mid_index
        rest = np.setdiff1d(np.arange(len(self.free)), [k])
        flat, _, _ = self.flat_state()
        out = np.empty((len(w_values), len(self.free)))
        lam = np.empty(len(w_values))
        for sign in (1, -1):
            idx = [i for i in np.argsort(np.abs(w_values)) if np.sign(w_values[i]) in (sign, 0)]
            q = flat.copy()
            prev = 0.0
            for i in idx:
                w = w_values[i]
                # small sub-steps for robust continuation
                for ws in np.linspace(prev, w, 4)[1:]:
                    q[k] = ws
                    for _ in range(maxiter):
                        f, K = self.force_tangent(q)
                        if np.linalg.norm(f[rest]) < tol*max(1.0, self.EA*self.end_shortening):
                            break
                        q[rest] -= np.linalg.solve(K[np.ix_(rest, rest)], f[rest])
                    else:
                        raise RuntimeError(f"constrained Newton failed at w={ws}")
                prev = w
                out[i] = q
                _, K = self.force_tangent(q)
                lam[i] = np.linalg.eigvalsh(K[np.ix_(rest, rest)])[0]
        return out, lam

    def classify(self, q, v, barrier):
        """+1 up, -1 down, 0 unresolved (total energy >= barrier). Batched."""
        e = self.energy(q, v)
        w = np.asarray(q)[..., self.mid_index]
        out = np.zeros(e.shape, dtype=np.int8)
        trapped = np.isfinite(e) & (e < barrier)
        out[trapped & (w > 0)] = 1
        out[trapped & (w < 0)] = -1
        return out
