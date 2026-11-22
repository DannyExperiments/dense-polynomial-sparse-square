#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 check_inputs.py "$@"
python3 - <<'PY'
from pathlib import Path
import re
paths = [*Path("SparseSquare").glob("*.lean"), Path("SparseSquare.lean"), Path("Audit.lean")]
for path in paths:
    source = path.read_text()
    if re.search(r"\b(sorry|admit|unsafe|native_decide)\b|^\s*(axiom|constant)\b", source, re.M):
        raise SystemExit("Prohibited Lean declaration or tactic in " + str(path))
PY
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT
cp -R SparseSquare SparseSquare.lean Audit.lean lakefile.toml lake-manifest.json lean-toolchain "$build_dir/"
cd "$build_dir"
lake env lean --version
lake build
lake env leanchecker SparseSquare.Endpoint SparseSquare.Certificate
lake env lean Audit.lean > audit.log
cat audit.log
python3 - <<'PY'
from pathlib import Path
import re
decls = re.findall(r"^#print axioms ([\w.]+)$", Path("Audit.lean").read_text(), re.M)
results = {}
for line in Path("audit.log").read_text().splitlines():
    empty = re.fullmatch("'([\\w.]+)' does not depend on any axioms", line)
    some = re.fullmatch("'([\\w.]+)' depends on axioms: \\[(.*)\\]", line)
    if empty:
        results[empty[1]] = set()
    elif some:
        results[some[1]] = set(some[2].split(", "))
    else:
        raise SystemExit("Unexpected audit output: " + line)
if set(results) != set(decls):
    raise SystemExit("An audited declaration was missing")
for name, axioms in results.items():
    allowed = set() if name.startswith("SparseSquare.Certificate.") else {"propext", "Quot.sound"}
    if not axioms <= allowed:
        raise SystemExit("Unexpected axioms for " + name)
print(f"Checked axiom dependencies for all {len(results)} declarations.")
PY
