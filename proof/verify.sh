#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
python3 certify.py
cmp proth-periods.csv checks/sparse/proth-periods.csv
python3 checks/construct/check_local_mask.py --tamper
python3 checks/arithmetic/sparse_repair_control.py
python3 checks/sparse/local_sections.py --test 7 13 23 37 101
python3 checks/sparse/fingerprint.py 1 2 17 314159 | diff - fingerprints.txt
python3 checks/sparse/verify_masks.py
python3 verify_compressed.py
