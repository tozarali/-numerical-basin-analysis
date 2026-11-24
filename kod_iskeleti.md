# kod_iskeleti.md — basin_project kod iskeletleri

> `opus.md`'nin eki. VS Code asistanı yeni modül eklerken bu dosyayı **sözleşme** olarak kullanır.
> Bu dosyadaki kodlar mevcut `basinlab` arayüzüne (`model.force`, `integrators.parameters`, `integrators.step`, `metrics.*`) uyumludur ve **24 Kasım 2026'da mevcut paketle birlikte çalıştırılıp test edildi: 28/28 test geçti** (13 eski + 15 yeni, ~25 sn).
> Durum etiketleri: **[TESTLİ]** çalışır ve testi var · **[İSKELET]** arayüz + kısmi uygulama, TODO'lar tamamlanmadan sonuç üretmez · **[HİPOTEZ]** çalışır ama bilimsel geçerliliği henüz sınanmadı.

---

## 0. Kullanım (asistan için)

1. Bu dosyadaki her modülü `basin_project/basinlab/` altına aynı adla yaz (veya ekteki zip'i aç). Testleri `basin_project/tests/test_skeletons.py` olarak ekle.
2. `python -m unittest discover -s tests -v` → hepsi geçmeli. Geçmeyen modülü kullanma.
3. **[İSKELET]** modüllerdeki TODO'ları sırayla tamamla; her TODO için önce test yaz.
4. **[HİPOTEZ]** modüllerin (D_sep, D_th) çıktıları yalnız `status="exploratory"` ile raporlanır; doğrulayıcı kullanım `protocol/PREREG.md` commit edildikten sonra.
5. Var olan dosyaları (integrators.py, metrics.py, model.py, reference.py, truss.py, experiments/*) **değiştirme**; yeni modüller onları içe aktarır.

## 1. Hedef klasör yapısı

```
paper2_idea/
├─ opus.md, kod_iskeleti.md, STATUS.md, DECISIONS.md
├─ IDEA.md, PROMPT.md, pilots/            (tarihsel, dokunma)
└─ basin_project/
   ├─ requirements.txt                    numpy, scipy, matplotlib (sürümler sabit)
   ├─ protocol/PREREG.md                  doğrulayıcı deneylerden önce commit
   ├─ basinlab/
   │  ├─ model.py integrators.py reference.py metrics.py truss.py   (mevcut)
   │  ├─ runmeta.py          [TESTLİ]   summary.json kaynak bilgisi
   │  ├─ linear_analysis.py  [TESTLİ]   amplifikasyon, ξ̄, periyot hatası
   │  ├─ dtheory.py          [HİPOTEZ]  T(E), D_sep, D_th, crossing
   │  ├─ stats.py            [TESTLİ]   parabolik ρ*, bootstrap GA, Wilson GA
   │  ├─ sampling.py         [TESTLİ]   Sobol IC, enerji bandı
   │  ├─ metrics_ext.py      [TESTLİ*]  net transfer, işaretli kayma, kutu taraması, belirsizlik üssü
   │  ├─ integrators_extra.py[TESTLİ/İSKELET] discrete gradient, BDF2, ρ∞-Bathe(TODO)
   │  ├─ beam2mode.py        [TESTLİ]   iki modlu burkulmuş kiriş
   │  ├─ mdof_genalpha.py    [TESTLİ]   çok DOF gen-α (batch)
   │  ├─ fe_arch.py          [İSKELET]  korotasyonel kiriş, sığ kemer
   │  ├─ diagnostic.py       [TESTLİ]   iki-koşu tanı testi
   │  ├─ experiments/_template.py        deney şablonu
   │  └─ figures/_template.py            figür şablonu
   ├─ tests/test_skeletons.py
   ├─ data/<deney>/summary.json + ham dosyalar
   ├─ docs/  manuscript/  figures/  zenodo/
```
`*` metrics_ext içindeki `signed_boundary_shift` ve `uncertainty_exponent` için henüz test yok: ilk iş sentetik bir sınır (düz çizgi kaydırılmış) ile test yaz.

## 2. Ortak sözleşmeler

| Konu | Kural |
|---|---|
| Etiketler | `int8`, {−1, 0, +1}; 0 = çözülemedi (asla üçüncü kuyu sayılmaz) |
| Integratör parametresi | `p = (am, af, gamma, beta)`; ailelere `integrators.parameters_for_spectral_radius(family, rho)` ile ulaşılır |
| Çok DOF modeli | `model.force_tangent(q) -> (f[...,n], K[...,n,n])`, `model.energy(q,v)`, `model.classify(q,v)`, `model.barrier` |
| Zaman | Her karşılaştırma **aynı bitiş zamanında**; son adım kısaltılır (`integrate_to`). BDF2 istisna: `end/h` tam sayı olmalı |
| Batch | IC'ler ön eksende; Newton yakınsaması tüm batch için kontrol edilir |
| Kayıt | Her deney `runmeta.experiment(name, config, status=...)` içinde; seed, commit, sürüm, süre otomatik |
| Hata ölçüsü | `errors[n_rho, n_ic]` bool matrisi sakla; tüm istatistikler (bootstrap dahil) bundan türetilir |
| Rastgelelik | Tüm seed'ler config'te |

**Not (ρ∞ eşlemesi):** gen-α'nın spektral yarıçapı Ωh=1e3'te 0.407 (nominal 0.4), Ωh→∞'da 0.400'e iner; HHT/WBZ/Newmark Ωh=1e3'te zaten limitte. Mevcut eşleme bu nedenle tutarlı, ama makalede "ρ∞ = Ωh→∞ limiti" diye tanımla ve `linear_analysis.spectral_radius(p, 1e8)` ile raporla.


## 3. runmeta.py [TESTLİ]

`basin_project/basinlab/runmeta.py`

Her deneyin `summary.json`'una tek tip kaynak bilgisi yazar; gövde hata verirse de dosyayı `outcome: failed` ile bırakır.

```python
"""Uniform provenance block for every experiment summary.json."""
import json
import platform
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import scipy


def git_commit(root):
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                             capture_output=True, text=True, timeout=5)
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                               capture_output=True, text=True, timeout=5)
        if out.returncode:
            return None
        return out.stdout.strip() + ("+dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return None


def environment(root):
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "scipy": scipy.__version__, "platform": platform.platform(),
            "git_commit": git_commit(root)}


@contextmanager
def experiment(name, config, *, status="exploratory"):
    """Usage:
        with experiment("dsep_exploratory", cfg) as run:
            ...; run["results"] = {...}
    Writes data/<name>/summary.json even if the body raises (status=failed).
    status must be 'exploratory' or 'confirmatory' (confirmatory requires a
    committed protocol/PREREG.md; see opus.md section 0)."""
    if status not in {"exploratory", "confirmatory"}:
        raise ValueError("status must be exploratory or confirmatory")
    root = Path(__file__).resolve().parents[1]
    out = root / "data" / name
    out.mkdir(parents=True, exist_ok=True)
    run = {"name": name, "status": status, "configuration": config,
           "environment": environment(root), "results": {}, "out_dir": str(out)}
    started = time.perf_counter()
    try:
        yield run
        run["outcome"] = "completed"
    except BaseException as error:
        run["outcome"] = f"failed: {type(error).__name__}: {error}"
        raise
    finally:
        run["runtime_seconds"] = round(time.perf_counter() - started, 3)
        (out / "summary.json").write_text(json.dumps(run, indent=2, default=float))
```

## 4. linear_analysis.py [TESTLİ]

`basin_project/basinlab/linear_analysis.py`

`integrators.step` ile aynı denklemlerden kurulan 3×3 amplifikasyon matrisi. WP3-1 ve 'iki bileşen' (sönüm + faz) hipotezinin temeli. Testler: aile ρ∞ eşlemesi, trapez kapalı formu (ξ̄=0, ΔT/T = Ωh/(2·atan(Ωh/2)) − 1), küçük adımda fiziksel ξ'nin geri gelmesi.

```python
"""Exact linear amplification analysis of the alpha/Newmark family.

Linear test problem: x'' + 2 xi W x' + W^2 x = 0, step h, Omega = W h.
State (x, v, a) -> (x1, v1, a1) through the same equations as
integrators.step, so results apply directly to `parameters(...)` tuples.
"""
import numpy as np


def amplification_matrix(p, omega_h, xi=0.0):
    """3x3 map for unit frequency (W=1) and step omega_h."""
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
    lam = np.linalg.eigvals(amplification_matrix(p, omega_h, xi))
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
```

## 5. dtheory.py [HİPOTEZ]

`basin_project/basinlab/dtheory.py`

Pencereden bağımsız D adayları. `d_sep`: eyer noktasından E₀=E_b(1+δ) ile başlayan **tek bir** kuyu gezintisinde (yarım koruyucu periyot) sayısal/referans enerji kaybı oranı. `d_theory`: seviye-0 ortalama teori (orbit harmonik varsayılır). Testler yalnızca h→0'da D→1 olduğunu doğrular; **ρ\* tahmin gücü sınanmadı.** Seviye-1 TODO: kayıpları orbit boyunca c·v² ile ağırlıklandır.

```python
"""Horizon-free energy-loss ratios: D_sep (one near-separatrix excursion,
measured) and D_th (averaged linear theory). HYPOTHESIS CODE: both are
candidates to be tested against observed basin optima (opus.md 3.2)."""
import numpy as np
from scipy.integrate import solve_ivp

from .integrators import step
from .model import energy, force
from .linear_analysis import principal_pair

BARRIER = 0.25  # double well V=(x^2-1)^2/4, saddle at x=0


def _gen_alpha_to(x, v, end, h, damping, p):
    """Scalar/array gen-alpha with an exact shortened last step."""
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    a = -damping*v - force(x)
    n = int(np.floor(end/h + 1e-12))
    for _ in range(n):
        x, v, a = step(x, v, a, h, damping, p)
    rem = end - n*h
    if rem > 1e-10:
        x, v, a = step(x, v, a, rem, damping, p)
    return x, v


def _reference_to(x, v, end, damping, rtol=1e-11, atol=1e-13):
    sol = solve_ivp(lambda t, y: (y[1], -damping*y[1] - force(y[0])),
                    (0, end), (x, v), method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0, -1], sol.y[1, -1]


def conservative_period(E):
    """Period of the undamped orbit with energy E (double well)."""
    if E <= 0 or abs(E - BARRIER) < 1e-9:
        raise ValueError("E must be positive and not on the separatrix")
    if E < BARRIER:
        x0, v0 = 1.0, np.sqrt(2*E)          # well bottom, moving right
    else:
        x0, v0 = 0.0, np.sqrt(2*(E - BARRIER))  # saddle, moving right
    ev = lambda t, y: y[0] - x0
    ev.direction = 1
    sol = solve_ivp(lambda t, y: (y[1], -force(y[0])), (0, 1e4), (x0 - 1e-12, v0),
                    method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    times = sol.t_events[0]
    times = times[times > 1e-6]
    if not len(times):
        raise RuntimeError("period not found")
    return float(times[0])


def d_sep(p, h, damping, delta=0.05, window="half_period"):
    """Measured energy-loss ratio over ONE excursion starting at the saddle
    with E0 = BARRIER*(1+delta). Window = half the conservative period of E0
    (one well excursion), matched exactly in time for both solvers."""
    E0 = BARRIER * (1 + delta)
    x0, v0 = 0.0, np.sqrt(2*(E0 - BARRIER))
    tau = conservative_period(E0) / 2 if window == "half_period" else float(window)
    xr, vr = _reference_to(x0, v0, tau, damping)
    xn, vn = _gen_alpha_to(x0, v0, tau, h, damping, p)
    d_ref = E0 - energy(xr, vr)
    d_num = E0 - energy(xn, vn)
    return {"D_sep": float(d_num/d_ref), "E0": E0, "window": tau,
            "dE_ref": float(d_ref), "dE_num": float(d_num)}


def d_theory(p, h, damping, E):
    """Level-0 averaged theory: treat the orbit of energy E as harmonic with
    frequency w(E)=2*pi/T(E); ratio of numerical to physical energy loss per
    cycle from the exact linear amplification at (w h, xi=c/(2w)).
    KNOWN LIMITATION: near the separatrix the orbit is far from harmonic;
    level-1 (TODO) weights local losses by c*v^2 along the orbit."""
    w = 2*np.pi / conservative_period(E)
    xi = damping / (2*w)
    lam = principal_pair(p, w*h, xi)
    if lam is None:
        return np.nan
    steps_per_cycle = 2*np.pi / (w*h)
    retained_num = abs(lam) ** (2*steps_per_cycle)
    retained_phys = np.exp(-2*xi*w * 2*np.pi/w)   # = exp(-2*pi*c/w)
    return float((1 - retained_num) / (1 - retained_phys))


def crossing(rhos, values, target=1.0):
    """Linear interpolation of the first crossing of `target` (None if none)."""
    rhos, values = np.asarray(rhos, float), np.asarray(values, float)
    idx = np.flatnonzero((values[:-1]-target)*(values[1:]-target) <= 0)
    if not len(idx):
        return None
    i = int(idx[0])
    if values[i+1] == values[i]:
        return float((rhos[i]+rhos[i+1])/2)
    return float(rhos[i] + (target-values[i])*(rhos[i+1]-rhos[i])/(values[i+1]-values[i]))
```

## 6. stats.py [TESTLİ]

`basin_project/basinlab/stats.py`

Optimum yerinin sürekli tahmini (parabolik tepe) ve IC'lerin ortak yeniden örneklemesiyle bootstrap %95 GA. Bayrak `ok` değilse (kenar, dışbükey değil) sonuç raporda işaretlenir.

```python
"""Continuous optimum location and bootstrap uncertainty over ICs."""
import numpy as np


def parabolic_minimum(rhos, values, half_width=3):
    """Fit a parabola to the points around the discrete minimum.
    Returns the vertex if it is a minimum inside the fitted window, else the
    discrete argmin (flagged)."""
    rhos, values = np.asarray(rhos, float), np.asarray(values, float)
    i = int(np.nanargmin(values))
    lo, hi = max(0, i-half_width), min(len(rhos), i+half_width+1)
    if hi - lo < 3:
        return float(rhos[i]), "edge"
    a, b, _ = np.polyfit(rhos[lo:hi], values[lo:hi], 2)
    if a <= 0:
        return float(rhos[i]), "not_convex"
    vertex = -b/(2*a)
    if not rhos[lo] <= vertex <= rhos[hi-1]:
        return float(rhos[i]), "outside"
    return float(vertex), "ok"


def bootstrap_optimum(rhos, errors, n_boot=1000, seed=0, half_width=3):
    """errors: bool array (n_rho, n_ic), True = label differs from reference
    for that IC. ICs are resampled jointly across rho (they are shared),
    so the CI reflects IC-sampling uncertainty of the optimum location."""
    errors = np.asarray(errors, bool)
    rng = np.random.default_rng(seed)
    n_ic = errors.shape[1]
    point, flag = parabolic_minimum(rhos, errors.mean(axis=1), half_width)
    samples = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n_ic, n_ic)
        samples[b] = parabolic_minimum(rhos, errors[:, idx].mean(axis=1), half_width)[0]
    lo, hi = np.percentile(samples, [2.5, 97.5])
    return {"rho_star": point, "flag": flag, "ci95": [float(lo), float(hi)],
            "boot_sd": float(samples.std(ddof=1))}


def proportion_ci(k, n, z=1.96):
    """Wilson interval for a mismatch fraction."""
    if n == 0:
        return (np.nan, np.nan)
    p = k/n
    den = 1 + z*z/n
    centre = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return (float(centre-half), float(centre+half))
```

## 7. sampling.py [TESTLİ]

`basin_project/basinlab/sampling.py`

PREREG'de sabitlenecek IC popülasyonları. Öneri: kutu Sobol (N≥2^15) + birincil popülasyon `energy_band(E_b, k·E_b]`.

```python
"""Pre-registrable IC populations (fixed seed, recorded in summary.json)."""
import numpy as np
from scipy.stats import qmc


def sobol_box(n, lower, upper, seed=12345):
    """Scrambled Sobol points in a box; n is rounded up to a power of 2."""
    m = int(np.ceil(np.log2(n)))
    pts = qmc.Sobol(d=len(lower), scramble=True, seed=seed).random_base2(m)
    return qmc.scale(pts, lower, upper)


def energy_band(points, energy_fn, e_low, e_high):
    """Keep points with e_low < E0 <= e_high (e.g. (E_b, k*E_b])."""
    e = energy_fn(*points.T)
    keep = (e > e_low) & (e <= e_high)
    return points[keep], e[keep]
```

## 8. metrics_ext.py [TESTLİ*]

`basin_project/basinlab/metrics_ext.py`

Doygunluğa dirençli ölçüler. `net_basin_transfer` işaretli alan ölçüsü (hangi havza büyüdü). `signed_boundary_shift` ve `uncertainty_exponent` için test eksik.

```python
"""Additional, saturation-resistant basin measures (opus.md WP1-3)."""
import numpy as np

from .metrics import boundary_points, entropy


def net_basin_transfer(reference, candidate):
    """Signed area measure: fraction of resolved ICs moved -1 -> +1 minus
    +1 -> -1. Zero for symmetric errors; sign shows which basin grows."""
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    valid = (reference != 0) & (candidate != 0)
    n = int(valid.sum())
    if n == 0:
        return {"net": None, "to_plus": None, "to_minus": None, "valid": 0}
    to_plus = np.mean((reference == -1) & (candidate == 1) & valid) * valid.size / n
    to_minus = np.mean((reference == 1) & (candidate == -1) & valid) * valid.size / n
    return {"net": float(to_plus - to_minus), "to_plus": float(to_plus),
            "to_minus": float(to_minus), "valid": n}


def signed_boundary_shift(reference, candidate, spacing):
    """For each reference boundary midpoint, the distance to the nearest
    candidate boundary midpoint, signed +1 if the candidate label at the
    nearer reference-side cell is +1 (the + basin grew), -1 otherwise.
    Returns the mean signed and mean absolute shift. Grid units = spacing."""
    from scipy.spatial import cKDTree
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    if np.any(reference == 0) or np.any(candidate == 0):
        raise ValueError("requires fully resolved grids")
    a = boundary_points(reference, spacing)
    b = boundary_points(candidate, spacing)
    if not len(a) or not len(b):
        return {"mean_signed": None, "mean_abs": None}
    dist, _ = cKDTree(b).query(a)
    dy, dx = spacing
    # sign: sample the candidate label at the reference boundary point's cell
    rows = np.clip(np.round(a[:, 0]/dy).astype(int), 0, reference.shape[0]-1)
    cols = np.clip(np.round(a[:, 1]/dx).astype(int), 0, reference.shape[1]-1)
    changed = candidate[rows, cols] != reference[rows, cols]
    sign = np.where(changed, candidate[rows, cols], 0)  # unchanged cell -> sign 0
    return {"mean_signed": float(np.mean(sign*dist)), "mean_abs": float(dist.mean()),
            "points": int(len(a))}


def entropy_box_sweep(labels, box_sizes=(2, 4, 5, 8, 10, 20)):
    """Sb/Sbb for several box widths (only those dividing the grid)."""
    n = np.asarray(labels).shape[0]
    return {int(b): entropy(labels, b) for b in box_sizes if n % b == 0}


def uncertainty_exponent(classify_ics, ics, eps_values, seed=0):
    """Grebogi-McDonald-Ott-Yorke: f(eps) = fraction of ICs whose outcome
    changes under a random +/-eps perturbation; f ~ eps^alpha.
    classify_ics: callable (N,2) -> labels in {-1,0,1}. Expensive: each eps
    needs two classifications of all ICs (reuse the base labels)."""
    rng = np.random.default_rng(seed)
    base = classify_ics(ics)
    out = []
    for eps in eps_values:
        direction = rng.normal(size=ics.shape)
        direction /= np.linalg.norm(direction, axis=1, keepdims=True)
        pert = classify_ics(ics + eps*direction)
        ok = (base != 0) & (pert != 0)
        out.append(float(np.mean(base[ok] != pert[ok])) if ok.any() else np.nan)
    f = np.asarray(out)
    good = f > 0
    alpha = (float(np.polyfit(np.log(np.asarray(eps_values)[good]), np.log(f[good]), 1)[0])
             if good.sum() >= 3 else None)
    return {"eps": list(map(float, eps_values)), "f": out, "alpha": alpha}
```

## 9. integrators_extra.py [TESTLİ / İSKELET]

`basin_project/basinlab/integrators_extra.py`

Enerji-tutarlı kontrol (discrete gradient; c=0'da enerji 1e−9 göreli hassasiyetle korunuyor) ve BDF2 (ρ∞=0, 2. derece doğrulandı). **ρ∞-Bathe bilerek boş:** katsayılar orijinal makaleden alınmalı, ezberden yazılmamalı; önce linear_analysis tarzı testle ρ∞ doğrulanmalı. BDF2 için ayrıca `linear_analysis` benzeri birinci-mertebe amplifikasyon analizi TODO.

```python
"""Extra SDOF schemes for the double well (vectorised over ICs).
Same signature idea as integrators.step: state arrays in, state arrays out.
Each scheme is a 'point' (no rho parameter) or has its own parameter."""
import numpy as np

from .model import force, tangent


def potential(x):
    return 0.25*(x*x - 1)**2


def discrete_gradient_step(x, v, h, damping, tol=1e-12, maxiter=30):
    """Gonzalez/Itoh-Abe discrete gradient (energy-consistent control):
    x1-x0 = h(v0+v1)/2,  v1-v0 = -h[c(v0+v1)/2 + (V(x1)-V(x0))/(x1-x0)].
    With c=0 energy is conserved to solver tolerance."""
    x1 = x + h*v
    for _ in range(maxiter):
        dx = x1 - x
        small = np.abs(dx) < 1e-9
        dg = np.where(small, force(0.5*(x+x1)),
                      (potential(x1) - potential(x)) / np.where(small, 1, dx))
        v1 = 2*dx/h - v
        r = v1 - v + h*(damping*0.5*(v + v1) + dg)
        if np.all(np.abs(r) < tol):
            return x1, v1
        # derivative of r wrt x1 (dg derivative by finite difference, robust)
        e = 1e-7
        x1e = x1 + e
        dxe = x1e - x
        dge = np.where(np.abs(dxe) < 1e-9, force(0.5*(x+x1e)),
                       (potential(x1e) - potential(x)) / np.where(np.abs(dxe) < 1e-9, 1, dxe))
        drdx = 2/h + h*(damping/h + (dge - dg)/e)
        x1 = x1 - r/drdx
    raise RuntimeError("discrete gradient Newton failed")


def bdf2_integrate(x, v, end, h, damping, tol=1e-12, maxiter=30):
    """BDF2 on (x, v) with one trapezoidal start step; rho_inf = 0, O(h^2).
    Exact end time: requires end/h integer (keeps BDF2 coefficients fixed)."""
    n = round(end/h)
    if not np.isclose(n*h, end, atol=1e-10):
        raise ValueError("end must be a multiple of h for BDF2")
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    # start: trapezoidal (Newmark 1/4,1/2) single step via discrete Newton
    xm, vm = x.copy(), v.copy()
    x1 = x + h*v
    for _ in range(maxiter):
        v1 = 2*(x1 - x)/h - v
        r = (v1 - v) + h*0.5*(damping*(v+v1) + force(x) + force(x1))
        if np.all(np.abs(r) < tol):
            break
        x1 -= r / (2/h + h*0.5*(damping*2/h + tangent(x1)))
    x, v = x1, 2*(x1 - xm)/h - vm
    for _ in range(n - 1):
        xg = 2*x - xm
        for _ in range(maxiter):
            v1 = (3*xg - 4*x + xm) / (2*h)
            r = 3*v1 - 4*v + vm + 2*h*(damping*v1 + force(xg))
            if np.all(np.abs(r) < tol):
                break
            drdx = 9/(2*h) + 2*h*(damping*3/(2*h) + tangent(xg))
            xg = xg - r/drdx
        else:
            raise RuntimeError("BDF2 Newton failed")
        xm, vm, x, v = x, v, xg, (3*xg - 4*x + xm)/(2*h)
    return x, v


def bathe_rho_integrate(x, v, end, h, damping, rho_inf):
    """rho_inf-Bathe composite scheme (Noh & Bathe 2019).
    TODO: implement from the original paper; verify rho_inf by
    linear_analysis-style amplification of a linear SDOF before use.
    Do NOT write coefficients from memory."""
    raise NotImplementedError("Implement from Noh & Bathe (2019) and test first")
```

## 10. beam2mode.py [TESTLİ]

`basin_project/basinlab/beam2mode.py`

İkinci bağımsız yapısal model (WP4). Kuvvet ve Hessian sonlu farkla, denge noktaları ve kararlılıkları test edildi. Bariyer düz durumdaki eyer (V=0). Bağlaşım gücü α ve P ile taranır; kafes için kullanılan bağlaşım diagnostiklerini burada da raporla.

```python
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
```

## 11. mdof_genalpha.py [TESTLİ]

`basin_project/basinlab/mdof_genalpha.py`

Kafes, kiriş ve FE için ortak çok DOF gen-α (IC batch'li). Beam2Mode'da DOP853'e karşı 2. derece yakınsama doğrulandı.

```python
"""Generalized-alpha for M a + C v + f(q) = 0, batched over ICs.
model.force_tangent(q) -> (f [...,n], K [...,n,n]); M, C constant (n,n).
Used for truss, Beam2Mode and the FE arch (reduce batch size for FE)."""
import numpy as np


def step(q, v, a, h, M, C, model, p, tol=1e-10, maxiter=30):
    am, af, gamma, beta = p
    c1 = 1/(beta*h*h)
    f0, _ = model.force_tangent(q)
    q1 = q + h*v + 0.5*h*h*a
    for it in range(maxiter + 1):
        a1 = c1*(q1 - q - h*v) - (0.5/beta - 1)*a
        v1 = v + h*((1-gamma)*a + gamma*a1)
        f1, K1 = model.force_tangent(q1)
        r = ((1-am)*a1 + am*a) @ M.T + ((1-af)*v1 + af*v) @ C.T + (1-af)*f1 + af*f0
        scale = max(1.0, float(np.max(np.abs(f0))))
        if np.max(np.abs(r)) <= tol*scale:
            return q1, v1, a1
        if it == maxiter:
            break
        J = (1-am)*c1*M + (1-af)*gamma*h*c1*C + (1-af)*K1
        q1 = q1 - np.linalg.solve(J, r[..., None])[..., 0]
    raise RuntimeError("mdof generalized-alpha Newton failed")


def integrate_to(q, v, end, h, M, C, model, p):
    """Exact end time via one shortened final step (as energy_loss_retest)."""
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    f, _ = model.force_tangent(q)
    a = np.linalg.solve(M, (-(v @ C.T) - f)[..., None])[..., 0]
    n = int(np.floor(end/h + 1e-12))
    for _ in range(n):
        q, v, a = step(q, v, a, h, M, C, model, p)
    rem = end - n*h
    if rem > 1e-10:
        q, v, a = step(q, v, a, rem, M, C, model, p)
    return q, v
```

## 12. fe_arch.py [İSKELET → GÜNCELLENDİ 24.11 öğleden sonra]

> **Not:** Aşağıdaki kod ilk iskelettir. Güncel `fe_arch.py` batch'li, analitik tanjantlı ve uç kısaltmalı burkulmuş kiriş modelidir (gerilmesiz sığ kemer yüksüz durumda bistabil çıkmadı). Güncel sürüm için dosyanın kendisine ve `CALISTIRMA.md`'ye bakın.

`basin_project/basinlab/fe_arch.py`

WP5 FE vakası. Eleman kinematiği (korotasyonel, Crisfield tipi) yazıldı; **rijit cisim hareketinde sıfır iç kuvvet** ve dinginlikte simetrik pozitif tanımlı tanjant test edildi. Eksikler (dosya başındaki 5 maddelik doğrulama listesi): analitik tanjant, statik denge/arc-length, eyer ve bariyer enerjisi, modal IC dilimi, `classify`, literatür benchmark karşılaştırması, openseespy çapraz kontrolü. Dönme ataleti katsayısı (`/78`) gerekçelendirilmeli veya tutarlı kütleye geçilmeli. Doğrulama bitmeden havza sonucu raporlanmaz.

```python
"""Shallow clamped arch with 2D corotational Euler-Bernoulli beam elements
(Crisfield-type). SKELETON: element kinematics implemented, tangent by
finite differences (replace with analytic tangent once tests pass).
DOFs per node: (u, w, theta). Clamped ends.

Validation required before any basin run (opus.md WP5):
  1. rigid-body rotation/translation gives zero internal force;
  2. small-deflection static response vs linear beam theory;
  3. static snap-through (limit) load vs a published shallow-arch benchmark
     (cite with DOI; do not quote values from memory);
  4. mesh refinement (n_el x2) of the barrier energy and of basin labels;
  5. a few dynamic trajectories cross-checked against openseespy.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Arch:
    span: float = 1.0
    rise: float = 0.05           # shallow: rise/span ~ 0.03-0.1
    n_el: int = 16
    EA: float = 1.0e3
    EI: float = 1.0
    rhoA: float = 1.0
    damping_mass: float = 0.0    # Rayleigh C = a0 M + a1 K0
    damping_stiff: float = 0.0

    def __post_init__(self):
        s = np.linspace(0, self.span, self.n_el + 1)
        self.X = np.column_stack((s, self.rise*np.sin(np.pi*s/self.span)))
        self.ndof = 3*(self.n_el + 1)
        fixed = [0, 1, 2, self.ndof-3, self.ndof-2, self.ndof-1]
        self.free = np.setdiff1d(np.arange(self.ndof), fixed)
        d = self.X[1:] - self.X[:-1]
        self.L0 = np.hypot(d[:, 0], d[:, 1])
        self.beta0 = np.arctan2(d[:, 1], d[:, 0])

    # --- element -----------------------------------------------------------
    def element_force(self, e, de):
        """de = (u1,w1,t1,u2,w2,t2) global displacements of element e."""
        X1, X2 = self.X[e], self.X[e+1]
        x1 = X1 + de[0:2]
        x2 = X2 + de[3:5]
        d = x2 - x1
        Ln = np.hypot(*d)
        beta = np.arctan2(d[1], d[0])
        rot = np.arctan2(np.sin(beta - self.beta0[e]), np.cos(beta - self.beta0[e]))
        t1, t2 = de[2] - rot, de[5] - rot          # local (deformational) rotations
        L0 = self.L0[e]
        N = self.EA*(Ln - L0)/L0
        M1 = self.EI/L0*(4*t1 + 2*t2)
        M2 = self.EI/L0*(2*t1 + 4*t2)
        c, s = np.cos(beta), np.sin(beta)
        B = np.array([[-c, -s, 0, c, s, 0],
                      [-s/Ln, c/Ln, 1, s/Ln, -c/Ln, 0],
                      [-s/Ln, c/Ln, 0, s/Ln, -c/Ln, 1]])
        return B.T @ np.array([N, M1, M2])

    def internal_force(self, d_full):
        f = np.zeros(self.ndof)
        for e in range(self.n_el):
            idx = np.arange(3*e, 3*e + 6)
            f[idx] += self.element_force(e, d_full[idx])
        return f

    def lumped_mass(self):
        m = np.zeros(self.ndof)
        for e, L in enumerate(self.L0):
            for node in (e, e+1):
                m[3*node:3*node+2] += 0.5*self.rhoA*L
                m[3*node+2] += self.rhoA*L**3/78     # rotary term: CHECK/justify
        return np.diag(m[self.free])

    # --- reduced (free-DOF) interface for mdof_genalpha ---------------------
    def force_tangent(self, q):
        """q: (..., n_free). Tangent by central differences (slow; TODO analytic)."""
        q = np.atleast_2d(q)
        n = len(self.free)
        F = np.empty_like(q)
        K = np.empty(q.shape + (n,))
        for b in range(q.shape[0]):
            full = np.zeros(self.ndof); full[self.free] = q[b]
            F[b] = self.internal_force(full)[self.free]
            for j in range(n):
                e = 1e-7*max(1.0, abs(q[b, j]))
                fp = full.copy(); fp[self.free[j]] += e
                fm = full.copy(); fm[self.free[j]] -= e
                K[b, :, j] = (self.internal_force(fp) - self.internal_force(fm))[self.free]/(2*e)
        return F, K

    # TODO: static equilibria (Newton + arc-length), saddle (minimum-mode
    # search) -> barrier energy; modal basis from K0, M for IC slices;
    # classify(q, v): energy < barrier and sign of midspan displacement.
```

## 13. diagnostic.py [TESTLİ]

`basin_project/basinlab/diagnostic.py`

WP6: referans olmadan iki ucuz koşunun (ρ=1 vs ρ=0.6 veya h vs h/2) uyuşmazlığı gerçek hatayı yakalıyor mu.

```python
"""Practical two-run diagnostic (opus.md WP6): without a reference, does the
label disagreement between two cheap runs predict the true error?"""
import numpy as np
from scipy.stats import spearmanr


def disagreement(labels_a, labels_b):
    a, b = np.asarray(labels_a), np.asarray(labels_b)
    ok = (a != 0) & (b != 0)
    return float(np.mean(a[ok] != b[ok])) if ok.any() else np.nan


def evaluate(cases, alarm=0.02, true_threshold=0.05):
    """cases: list of dicts with 'true_error' (vs reference) and
    'disagreement' (e.g. rho=1 vs rho=0.6, or h vs h/2).
    Returns rank correlation and alarm hit/false-alarm rates."""
    t = np.array([c["true_error"] for c in cases])
    d = np.array([c["disagreement"] for c in cases])
    big = t > true_threshold
    raised = d > alarm
    rho, pval = spearmanr(t, d)
    return {"spearman": float(rho), "n": int(len(t)),
            "hit_rate": float(raised[big].mean()) if big.any() else None,
            "false_alarm_rate": float(raised[~big].mean()) if (~big).any() else None,
            "alarm": alarm, "true_threshold": true_threshold,
            "note": "cases share models/ICs; p-values are descriptive only"}
```

## 14. tests/test_skeletons.py

`basin_project/tests/test_skeletons.py`

15 test; mevcut 13 testle birlikte ~25 sn.

```python
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from basinlab import (beam2mode, diagnostic, dtheory, fe_arch, integrators,
                      integrators_extra, linear_analysis, mdof_genalpha, stats)
from basinlab.model import energy, force


class LinearAnalysis(unittest.TestCase):
    def test_rho_infinity_matches_family_map(self):
        for fam, target in [("generalized-alpha", .4), ("HHT", .7), ("WBZ", .6), ("Newmark", .9)]:
            p, rho = integrators.parameters_for_spectral_radius(fam, target)
            self.assertAlmostEqual(linear_analysis.spectral_radius(p, 1e8), rho, places=3)

    def test_trapezoid_closed_form(self):
        p = integrators.parameters("generalized-alpha", 1.0)
        for w in (0.1, 0.5, 1.0):
            r = linear_analysis.dissipation_and_period(p, w)
            self.assertAlmostEqual(r["xi_bar"], 0.0, places=10)
            self.assertAlmostEqual(r["period_error"], w/(2*np.arctan(w/2)) - 1, places=10)

    def test_physical_damping_recovered_small_step(self):
        p = integrators.parameters("generalized-alpha", .5)
        r = linear_analysis.dissipation_and_period(p, 1e-3, xi=0.05)
        self.assertAlmostEqual(r["xi_bar"], 0.05, places=5)


class DTheory(unittest.TestCase):
    def test_small_energy_period(self):
        self.assertAlmostEqual(dtheory.conservative_period(1e-6), 2*np.pi/np.sqrt(2), places=4)

    def test_dsep_tends_to_one(self):
        p = integrators.parameters("generalized-alpha", .5)
        self.assertAlmostEqual(dtheory.d_sep(p, 0.01, 0.15)["D_sep"], 1.0, places=3)

    def test_dtheory_tends_to_one(self):
        p = integrators.parameters("generalized-alpha", .5)
        self.assertAlmostEqual(dtheory.d_theory(p, 1e-3, 0.15, 0.5), 1.0, places=3)

    def test_crossing(self):
        self.assertAlmostEqual(dtheory.crossing([0, 1], [2, 0]), 0.5)


class Stats(unittest.TestCase):
    def test_parabola_and_bootstrap(self):
        r = np.linspace(0, 1, 21)
        rng = np.random.default_rng(1)
        prob = 0.1 + 0.4*(r - 0.43)**2
        errors = rng.random((21, 20000)) < prob[:, None]
        out = stats.bootstrap_optimum(r, errors, n_boot=200)
        self.assertLess(abs(out["rho_star"] - 0.43), 0.05)
        self.assertLess(out["ci95"][0], out["rho_star"])


class ExtraIntegrators(unittest.TestCase):
    def test_discrete_gradient_conserves_energy(self):
        x, v = np.array([0.3, 1.5]), np.array([1.0, 0.2])
        e0 = energy(x, v)
        for _ in range(500):
            x, v = integrators_extra.discrete_gradient_step(x, v, 0.2, 0.0)
        np.testing.assert_allclose(energy(x, v), e0, rtol=1e-9)

    def test_bdf2_second_order(self):
        x0, v0, c, T = 0.5, 0.3, 0.15, 2.0
        ref = solve_ivp(lambda t, y: (y[1], -c*y[1] - force(y[0])), (0, T), (x0, v0),
                        method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1]
        errs = [abs(integrators_extra.bdf2_integrate(x0, v0, T, h, c)[0] - ref[0])
                for h in (0.02, 0.01)]
        self.assertAlmostEqual(np.log2(errs[0]/errs[1]), 2.0, delta=0.3)


class Beam(unittest.TestCase):
    def test_force_tangent_fd_and_equilibria(self):
        m = beam2mode.Beam2Mode(P=2.0, alpha=1.0)
        q = np.array([0.3, -0.2])
        f, K = m.force_tangent(q)
        e = 1e-6
        for j in range(2):
            dq = np.zeros(2); dq[j] = e
            np.testing.assert_allclose((m.potential(q+dq) - m.potential(q-dq))/(2*e), f[j], rtol=1e-6)
            np.testing.assert_allclose((m.force_tangent(q+dq)[0] - m.force_tangent(q-dq)[0])/(2*e),
                                       K[:, j], rtol=1e-6, atol=1e-8)
        for qe in m.equilibria():
            np.testing.assert_allclose(m.force_tangent(qe)[0], 0, atol=1e-12)
            self.assertTrue(np.all(np.linalg.eigvalsh(m.force_tangent(qe)[1]) > 0))

    def test_mdof_genalpha_converges(self):
        m = beam2mode.Beam2Mode()
        M, C = np.eye(2), np.diag(m.damping)
        q0, v0, T = np.array([[0.4, 0.1]]), np.array([[0.0, 0.3]]), 3.0
        def rhs(t, y):
            f, _ = m.force_tangent(y[:2]); return np.r_[y[2:], -C @ y[2:] - f]
        ref = solve_ivp(rhs, (0, T), np.r_[q0[0], v0[0]], method="DOP853",
                        rtol=1e-12, atol=1e-14).y[:2, -1]
        p = integrators.parameters("generalized-alpha", .6)
        errs = [np.linalg.norm(mdof_genalpha.integrate_to(q0, v0, T, h, M, C, m, p)[0][0] - ref)
                for h in (0.02, 0.01)]
        self.assertAlmostEqual(np.log2(errs[0]/errs[1]), 2.0, delta=0.3)


class Arch(unittest.TestCase):
    def test_rigid_body_motion_is_force_free(self):
        a = fe_arch.Arch(n_el=8)
        nodes = a.X
        th = 0.3
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        moved = nodes @ R.T + np.array([0.2, -0.1])
        d = np.zeros(a.ndof)
        d[0::3], d[1::3], d[2::3] = (moved - nodes)[:, 0], (moved - nodes)[:, 1], th
        np.testing.assert_allclose(a.internal_force(d), 0, atol=1e-9)

    def test_tangent_symmetric_at_rest(self):
        a = fe_arch.Arch(n_el=6)
        f, K = a.force_tangent(np.zeros(len(a.free)))
        np.testing.assert_allclose(f, 0, atol=1e-12)
        np.testing.assert_allclose(K[0], K[0].T, rtol=1e-5, atol=1e-3)
        self.assertTrue(np.all(np.linalg.eigvalsh(K[0]) > 0))


class Diagnostic(unittest.TestCase):
    def test_evaluate(self):
        cases = [{"true_error": t, "disagreement": t/2} for t in (0.01, 0.03, 0.08, 0.2)]
        out = diagnostic.evaluate(cases)
        self.assertAlmostEqual(out["spearman"], 1.0)
        self.assertEqual(out["hit_rate"], 1.0)


if __name__ == "__main__":
    unittest.main()
```
## 15. Deney şablonu — `basin_project/basinlab/experiments/_template.py` [TESTLİ]

Smoke testi (723 IC, 7 ρ, ~1 dk) uçtan uca çalıştı: Sobol örnekleme → enerji bandı → referans → aile taraması → `errors.npz` → bootstrap ρ\* → `summary.json`. Yeni her deney bu kalıbın kopyasıdır.

```python
"""TEMPLATE — copy to <name>.py. Run: python -m basinlab.experiments.<name>
Exploratory by default. Set STATUS='confirmatory' only if protocol/PREREG.md
is committed and this config matches it exactly."""
import numpy as np

from ..integrators import parameters_for_spectral_radius
from ..model import classify, energy
from ..reference import rk4
from ..runmeta import experiment
from ..sampling import energy_band, sobol_box
from ..stats import bootstrap_optimum
from .energy_loss_retest import exact_step_grid

NAME = "template_smoke"
STATUS = "exploratory"
CONFIG = {"family": "generalized-alpha", "damping": 0.15, "h": 0.4,
          "rhos": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], "T_basin": 60.0,
          "reference_h": 0.01, "n_sobol": 1024, "seed": 12345,
          "box": [[-2, -2], [2, 2]], "energy_band": [0.25, 2.0]}


def main(cfg=CONFIG):
    with experiment(NAME, cfg, status=STATUS) as run:
        pts = sobol_box(cfg["n_sobol"], *cfg["box"], seed=cfg["seed"])
        pts, _ = energy_band(pts, energy, *cfg["energy_band"])
        x, v = pts[:, 0], pts[:, 1]
        ref = classify(*rk4(x, v, cfg["T_basin"], cfg["reference_h"], cfg["damping"]))
        rhos = np.asarray(cfg["rhos"])
        errors = np.zeros((len(rhos), len(x)), bool)
        unresolved = []
        for i, rho in enumerate(rhos):
            p, _ = parameters_for_spectral_radius(cfg["family"], float(rho))
            lab = classify(*exact_step_grid(x, v, cfg["T_basin"], cfg["h"], cfg["damping"], p))
            ok = (ref != 0) & (lab != 0)
            errors[i] = ok & (ref != lab)
            unresolved.append(int((~ok).sum()))
        np.savez_compressed(f"{run['out_dir']}/errors.npz", rhos=rhos, errors=errors, ics=pts)
        run["results"] = {"n_ic": int(len(x)), "unresolved_per_rho": unresolved,
                          "mismatch": errors.mean(axis=1).tolist(),
                          "optimum": bootstrap_optimum(rhos, errors, n_boot=300)}
    return run


if __name__ == "__main__":
    print(main()["results"])
```

**Dikkat (smoke testinde görüldü):** U eğrisi asimetrik ve sığ; Δρ=0.1 ızgarada `half_width=3` parabolü tepeyi ayrık minimumdan (0.5) belirgin şekilde uzağa (0.57) çekti. Doğrulayıcı deneylerde optimum çevresinde Δρ=0.025 kullan, pencereyi dar tut (±0.075) ve ayrık argmin'i de raporla. Pencere genişliği PREREG'de sabitlenir.

## 16. Figür şablonu — `basin_project/basinlab/figures/_template.py` [TESTLİ]

Kural: figür script'i hesap yapmaz, yalnızca `data/` okur; vektörel PDF yazar; renk körü dostu palet.

```python
"""TEMPLATE — one script per figure: python -m basinlab.figures.<name>
Reads only data/<experiment>/ files (never recomputes), writes vector PDF
to figures/. Colour-blind-safe palette, sizes for a single journal column."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
COL1, COL2 = 3.35, 6.9          # inches; check the target journal's guide
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]


def main(experiment="template_smoke"):
    d = np.load(ROOT / "data" / experiment / "errors.npz")
    s = json.loads((ROOT / "data" / experiment / "summary.json").read_text())
    fig, ax = plt.subplots(figsize=(COL1, 2.4), constrained_layout=True)
    ax.plot(d["rhos"], 100*d["errors"].mean(axis=1), "o-", color=PALETTE[0])
    lo, hi = s["results"]["optimum"]["ci95"]
    ax.axvspan(lo, hi, color=PALETTE[1], alpha=.25, lw=0, label="95% CI of $\\rho^*$")
    ax.set_xlabel(r"$\rho_\infty$"); ax.set_ylabel("basin mismatch (%)")
    ax.legend(frameon=False, fontsize=7)
    out = ROOT / "figures" / f"{experiment}.pdf"
    fig.savefig(out); print(out)


if __name__ == "__main__":
    main()
```

## 17. PREREG şablonu — `basin_project/protocol/PREREG_TEMPLATE.md`

```markdown
# PREREG — <experiment name>   (commit BEFORE running; then do not edit)

- Date / commit of code used:
- Question & hypothesis (H#):
- Models and parameters: family list, rho grid (fine step near expected optimum), h values, c values, epsilon
- IC population: generator (sobol_box seed=…, N=…), box, energy band (E_b, k*E_b], primary vs secondary population
- Reference: method, step/tolerance, T_basin; check (T vs 2T, h vs h/2) on the same ICs
- Primary measure: mismatch on primary population. Secondary: Sb/Sbb (box sizes …), net transfer, signed shift
- Saturation rule: cells with mismatch >= 30 % OR Sb >= 0.9 ln2 are flagged and excluded from correlations
- Optimum estimator: stats.bootstrap_optimum(half_width=…, n_boot=…); report discrete argmin as well
- D definition: (T_D median | D_sep delta=… | D_th level …); crossing by dtheory.crossing
- Success criterion (gate): e.g. |rho*_pred - rho*_obs| <= 0.1 in >= 80 % of reliable cells AND better than nominal-rho baseline
- What would count as failure, and what we will write then:
- Expected runtime / memory:
```

## 18. STATUS.md ve DECISIONS.md başlangıç içeriği

```markdown
# STATUS (son güncelleme: YYYY-MM-DD)
| Faz | Durum | Kanıt dosyası | Not |
|---|---|---|---|
| 0 Literatür | açık | docs/literature.md | |
| 1 Altyapı | açık | data/phase1_results.json | 200-400 ızgara yok |
| 2 Rejim | kısmi | | |
| 3 D mekanizması | açık | data/family_optimum_refinement | pencere/popülasyon bağımlılığı |
| 4 Yapısal/FE | açık | data/truss_basins | kiriş + FE yok |
## Son oturumda yapılanlar
## Açık sorular (kullanıcıya)
## Sıradaki adım
```

```markdown
# DECISIONS
| Tarih | Karar | Gerekçe | Geri alınabilir mi |
|---|---|---|---|
```

## 19. İlk deney sırası (bu iskeletlerle)

1. **`linear_curves`** (saniyeler): 4 aile × ρ∞ ∈ {0.2…1} × Ωh ∈ [0.01, 3] için ξ̄ ve ΔT/T. Figür: ξ̄(Ωh) ve ΔT/T(Ωh), Newmark'ın O(h) sönümü görünür olmalı.
2. **`two_component_exploratory`** (dakikalar): mevcut `data/family_optimum_refinement/rho_refinement.csv` + linear_analysis → mismatch ~ |D−1| + |ΔT/T| regresyonu. Yeni koşu yok, sadece mevcut veri.
3. **`dsep_exploratory`** (≈10 dk): 4 aile, c=0.15, h ∈ {0.4, 0.5, 0.6}, δ ∈ {0.02, 0.05, 0.1}; D_sep=1 kesişimi vs mevcut optimumlar. İlk hızlı bakış (bu dosya hazırlanırken, tek nokta): gen-α, h=0.4, δ=0.05'te D_sep=1 kesişimi ρ≈0.54; aynı alt kümedeki gözlenen optimum 0.475, T_D=4 medyanı 0.43. Yani **D_sep tek başına mevcut tanımdan daha iyi değil**; δ duyarlılığı ve faz bileşeniyle birlikte değerlendirilmeli. Seviye-0 D_th bu noktada 1'i hiç kesmedi → seviye-1 (c·v² ağırlıklı) gerekli.
4. Bu üçüne göre `protocol/PREREG.md` → kullanıcı onayı → WP1 ızgara yakınsaması (uzun koşu; önce süre tahmini).
5. Paralelde: `metrics_ext` eksik testleri, `fe_arch` doğrulama listesi madde 1–3.

**Hız notu:** `parameters_for_spectral_radius` HHT/WBZ/Newmark için her çağrıda 60 adımlık bisection yapıyor; taramalarda sonucu önbelleğe al (`functools.lru_cache` ile sarmala veya ρ→parametre tablosunu bir kez hesapla).
