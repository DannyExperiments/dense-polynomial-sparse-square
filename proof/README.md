# A dense polynomial whose square has fewer than the square root of its degree in terms

The five-section manuscript specifies an integer polynomial of degree n, proves that all n+1 coefficients are nonzero, and proves the stronger estimate wt(P²) < (501/1000)√n. The degree has 5,798 decimal digits. The exact degree and a conservative integer upper bound for wt(P²) appear in certificate.json. No full coefficient array has been generated or submitted to Epoch. External peer review is outstanding.

## Read or replay

Open **sparse-squares.pdf** for the proof. **sparse-squares.tex** is its source. From this directory run the following with Python 3.9 or newer and a POSIX shell. No downloads or third-party Python packages are required.

    sh verify.sh

The script exits on failure. Use an ordinary Python invocation with assertions enabled. The exact verifier for the actual sparse factors can also run independently.

    python3 checks/sparse/verify_masks.py

For archive integrity use either of the following commands in this directory.

    shasum -a 256 -c SHA256SUMS
    sha256sum -c SHA256SUMS

## Contents and evidentiary scope

- **checks/sparse/compact-local-masks.json** contains the 128 local factors as signed sparse integer numerators. Each local polynomial is the recorded numerator divided by (1−z)^4. These divisions are exact. The manuscript combines those polynomials into a single dense witness. This data format is distinct from Epoch’s expanded-array submission format.
- **checks/sparse/verify_masks.py** reads the recorded numbers independently of their generator. On intervals between consecutive exponents of a squared numerator the coefficients of its quotient by (1−z)^8 are polynomial of degree at most seven. Checking eight points, or checking every point in a shorter interval, verifies every coefficient in that interval. This proves that whole forbidden residue classes vanish even at the actual prime periods. The program checks all local heights, all cyclic values and derivatives, and the final integer inequalities.
- **proth-periods.csv**, **proth.py**, **certify.py**, and **certificate.json** fix the periods and check their primality certificates, the deterministic selection rule, and the integer bounds in the text. The duplicate CSV alongside the local tools makes those tools directly runnable from any directory.
- **checks/sparse/local_sections.py** builds a section of a local factor from the rank-two system. **checks/sparse/export_masks.py** regenerates the full sparse file. Running that optional script writes its output beside the script. **fingerprints.txt**, **checks/sparse/fingerprint.py**, and **verify_compressed.py** provide two routes to evaluating the actual giant witness at the same four finite-field points, one from the section equations and one from the recorded numerators. The modular evaluations are reproducibility checks.
- **checks/construct/check_local_mask.py** independently builds and directly squares local factors at eleven small primes through 257. Its tamper control intentionally changes one coefficient and verifies that the checks catch it.
- **checks/arithmetic/sparse_repair_control.py**, **checks/director/two_mask.py**, and **checks/sparse/hermite.py** exercise the global density repair in a smaller two-prime example. Direct integer multiplication agrees with the decomposed square at degree 1,168. The small control has a dense square and does not meet the strict threshold. The 128-prime example is the one certified to cross it.
- **verification.log** records the outputs from a replay. **SHA256SUMS** lists the hashes of all other files in the archive.

## Sources

The [Epoch problem statement](https://epoch.ai/frontiermath/open-problems/polynomial-square-few-coefficients) specifies the degree convention and asks for an expanded array of integer coefficients. Mahmoud Elgamal’s [November 2026 preprint](https://doi.org/10.5281/zenodo.24128991) gives a related construction with asymptotic leading constant 4√3. Section 4 derives the auxiliary polynomial identities directly.
