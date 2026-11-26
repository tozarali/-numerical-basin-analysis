"""Exploratory: do two cheap runs without a reference flag real error?
Uses regime_map part files. Predictors:
  two-rho : disagreement between labels at rho and at rho=1 (same family)
  two-h   : disagreement between labels at h and at h/2 (same family, rho)
Target: true mismatch vs the RK4 reference. Saturated cells excluded."""
import json
from pathlib import Path

import numpy as np

from ..diagnostic import disagreement, evaluate
from ..runmeta import experiment
from ._common import ROOT, out_name, quick, saturated


def main():
    src = ROOT / "data" / out_name("regime_map") / "parts"
    if not src.exists():
        raise FileNotFoundError(f"run regime_map first ({src})")
    cfg = {"source": str(src.relative_to(ROOT)), "alarm": 0.02, "true_threshold": 0.05}
    with experiment(out_name("diagnostic_eval"), cfg) as run:
        data = {}
        for meta_file in sorted(src.glob("c*_h*.json")):
            meta = json.loads(meta_file.read_text(encoding="utf-8"))
            labels = np.load(meta_file.with_suffix(".npz"))["labels"]
            c, h = meta["cells"][0]["c"], meta["cells"][0]["h"]
            data[(c, h)] = (meta, labels)
        two_rho, two_h = [], []
        for (c, h), (meta, labels) in data.items():
            cases = [tuple(x) for x in meta["cases"]]
            index = {cs: i for i, cs in enumerate(cases)}
            for i, row in enumerate(meta["cells"]):
                fam, rho = row["family"], row["rho"]
                if rho is None or row["mismatch"] is None or row["saturated"]:
                    continue
                j = index.get((fam, 1.0))
                if j is not None and rho != 1.0:
                    two_rho.append({"c": c, "h": h, "family": fam, "rho": rho,
                                    "true_error": row["mismatch"],
                                    "disagreement": disagreement(labels[i], labels[j])})
                half = data.get((c, round(h/2, 6)))
                if half is not None:
                    hcases = {tuple(x): k for k, x in enumerate(half[0]["cases"])}
                    k = hcases.get((fam, rho))
                    if k is not None:
                        two_h.append({"c": c, "h": h, "family": fam, "rho": rho,
                                      "true_error": row["mismatch"],
                                      "disagreement": disagreement(labels[i], half[1][k])})
        res = {}
        for name, cases in (("two_rho", two_rho), ("two_h", two_h)):
            res[name] = evaluate(cases, cfg["alarm"], cfg["true_threshold"]) if len(cases) >= 3 \
                else {"n": len(cases), "note": "too few cases"}
        (Path(run["out_dir"]) / "cases.json").write_text(
            json.dumps({"two_rho": two_rho, "two_h": two_h}, indent=1), encoding="utf-8")
        run["results"] = res
    return run


if __name__ == "__main__":
    main()
