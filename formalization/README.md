# Lean verification

The project pins Lean 4.30.0 and uses its bundled standard library. Lake has no external dependency to resolve. With this toolchain and Python 3 installed, run from this directory:

```sh
./verify.sh
```

The script compares the Lean input literals with the JSON certificate and period table in `../proof/`, builds a fresh copy in a temporary directory, runs the bundled Lean kernel checker, and audits the axioms of the named theorems. The script accepts an alternate input directory through `--proof-dir PATH`.

## Formal statements

- `SparseSquare.middleThird_sumfree` proves that the open middle third is sum-free modulo every natural modulus. The proof splits according to whether adding two representatives wraps around the period.
- `SparseSquare.coefficients_zero_of_evaluate_zero` proves that a polynomial with signed integer coefficients of absolute value less than the absolute value of a nonzero integer argument can vanish there only if every coefficient is zero. The induction uses divisibility of the constant coefficient by the argument.
- `SparseSquare.scaledProductCoeff_ne_zero` proves the endpoint repair's density lemma. If the first and last coefficients of E are nonzero, V is dense through degree t, and t is at least the degree d assigned to E, all coefficients of the ordinary convolution E(x)V(zx) through d+t are nonzero when the integer scaling parameter exceeds the magnitudes of the individual products. The convolution is defined in coefficient form by Horner evaluation of each diagonal.
- `SparseSquare.Certificate` recomputes the product of the 128 listed periods, the auxiliary ceiling square root, the repetition count, the quotient support bound, endpoint separation, the 5,798-digit degree and the support bound of the paper. It verifies both `(1000 * weightUpper)^2 < 501^2 * degree` and `weightUpper^2 < degree`. These closed integer propositions reduce in Lean's kernel. The theorem `transfer_to_weight` applies both inequalities to any weight proved at most the recorded upper bound.

## Scope

The remaining mathematical connection between the concrete polynomial and the global support bound and exact degree is established in the accompanying paper. In particular the local rational-function construction, the primality and local polynomial certificates, the CRT count, and the full assembly of the witness are outside this Lean development. The certificate inequalities alone do not establish existence of a sparse square. The numerical theorems and transfer theorem have empty axiom dependencies. The general coefficient lemmas use the standard Lean principles `propext` and `Quot.sound`.
