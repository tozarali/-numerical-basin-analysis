"""Approve protocol/PREREG.md and commit it. Run from anywhere:
    python protocol/approve_prereg.py
Asks for an explicit 'evet' before changing anything. Creates the Git repo
if missing, sets the status line, commits the whole project (data and
*.npz are ignored by .gitignore). Then run: run_all.cmd --confirmatory
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # basin_project
REPO = ROOT.parent
PREREG = ROOT / "protocol" / "PREREG.md"


def git(*a, check=True):
    r = subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True)
    if check and r.returncode != 0:
        sys.exit(f"git {' '.join(a)} failed:\n{r.stdout}\n{r.stderr}")
    return r


def main():
    text = PREREG.read_text(encoding="utf-8")
    if "Status: READY FOR APPROVAL" not in text and "Status: APPROVED" not in text:
        sys.exit("PREREG.md has an unexpected status line; aborting.")
    print(text[:1500], "\n...\n")
    print("Bu protokolu onayliyor musunuz? Onaydan sonra esikler degistirilemez.")
    if input("Onaylamak icin 'evet' yazin: ").strip().lower() != "evet":
        sys.exit("Onaylanmadi, hicbir sey degismedi.")
    if not (REPO / ".git").exists():
        git("init")
    if git("config", "user.name", check=False).stdout.strip() == "":
        git("config", "user.name", "Muhammed Koleoglu")
    if git("config", "user.email", check=False).stdout.strip() == "":
        git("config", "user.email", "mkoleoglu@mku.edu.tr")
    PREREG.write_text(text.replace("Status: READY FOR APPROVAL", "Status: APPROVED"), encoding="utf-8")
    git("add", "-A")
    git("commit", "-m", "Approve preregistration v1.0 and freeze code")
    print("Tamam. Commit:", git("rev-parse", "HEAD").stdout.strip())


if __name__ == "__main__":
    main()
