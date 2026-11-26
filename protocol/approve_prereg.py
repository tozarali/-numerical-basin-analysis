"""Approve protocol/PREREG.md and commit it.

    python protocol/approve_prereg.py

Asks for an explicit 'yes' before changing anything. Initializes a Git
repository if none exists, sets the status line to APPROVED and commits the
project. The confirmatory experiment refuses to run unless the protocol file
is approved, tracked and unmodified. Git user.name / user.email must already
be configured.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "protocol" / "PREREG.md"


def git(*a, check=True):
    r = subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True)
    if check and r.returncode != 0:
        sys.exit(f"git {' '.join(a)} failed:\n{r.stdout}\n{r.stderr}")
    return r


def main():
    text = PREREG.read_text(encoding="utf-8")
    if "Status: READY FOR APPROVAL" not in text and "Status: APPROVED" not in text:
        sys.exit("PREREG.md has an unexpected status line; aborting.")
    print(text[:1500], "\n...\n")
    print("Approve this protocol? After approval the thresholds must not be changed.")
    if input("Type 'yes' to approve: ").strip().lower() != "yes":
        sys.exit("Not approved; nothing changed.")
    if not (ROOT / ".git").exists():
        git("init")
    if not git("config", "user.email", check=False).stdout.strip():
        sys.exit("Configure git user.name and user.email first.")
    PREREG.write_text(text.replace("Status: READY FOR APPROVAL", "Status: APPROVED"), encoding="utf-8")
    git("add", "-A")
    git("commit", "-m", "Approve preregistration and freeze code")
    print("Done. Commit:", git("rev-parse", "HEAD").stdout.strip())


if __name__ == "__main__":
    main()
