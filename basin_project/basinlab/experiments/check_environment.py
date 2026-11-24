"""Verify the preparation environment against the project requirements."""
import importlib
from pathlib import Path
import sys

from packaging.requirements import Requirement


def main():
    print("Python:", sys.executable, sys.version.split()[0], flush=True)
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("Python 3.12 is required")
    requirements = Path(__file__).resolve().parents[2] / "requirements.txt"
    for line in requirements.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        requirement = Requirement(line)
        module = importlib.import_module(requirement.name)
        version = module.__version__
        print(requirement.name, version, flush=True)
        if version not in requirement.specifier:
            raise RuntimeError(f"{requirement.name} {version} does not meet {requirement.specifier}")


if __name__ == "__main__":
    main()
