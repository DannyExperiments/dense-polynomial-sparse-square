# A dense integer polynomial with a sparse square

[![Repository verification](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/verify.yml)
[![Partial Lean verification](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/lean.yml/badge.svg?branch=main)](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/lean.yml)
[![PDF build](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/pdf.yml/badge.svg?branch=main)](https://github.com/DannyExperiments/dense-polynomial-sparse-square/actions/workflows/pdf.yml)

Can an integer polynomial of degree n have a nonzero coefficient in every position while its square has fewer than √n nonzero coefficients? **Yes.** The paper specifies one finite polynomial by 128 compact factors. Its square has strictly fewer than 0.501√n nonzero coefficients. The degree has 5,798 decimal digits.

[Read the paper (PDF)](paper/manuscript.pdf) · [LaTeX source](paper/manuscript.tex) · [Original question](evidence/epoch-problem.md) · [Epoch's problem and submission instructions](https://epoch.ai/frontiermath/open-problems/polynomial-square-few-coefficients)

## How the construction works

Each local factor rules out an interval of residue classes for the exponents in its ordinary square. Combining 128 coprime periods by the Chinese remainder theorem makes the surviving fraction small enough. Repeating the resulting coefficient pattern fills the interior of the final polynomial. An explicit correction fills both ends while keeping its contribution to the square within the support bound. The paper proves the bound using integer arithmetic. The exact number of surviving terms in the square is unknown.

## What you can check

| Material | Purpose |
| --- | --- |
| [Original audited proof and certificate](proof/) | The manuscript reviewed by the referee, the 128 sparse numerators, the exact degree and the supplied verification programs. |
| [Independent audit and source archives](evidence/) | The AI referee's report, separately written checker and recorded results. Its mathematical verdict concerns the archived manuscript. |
| [Lean formalization](formalization/) | General coefficient and residue lemmas, and exact arithmetic for the final certificate. The paper supplies the derivation of the support bound for the constructed polynomial. |
| [Interactive illustration](visualization/) | Download the HTML file and open it locally to explore the residue restrictions and endpoint correction. |

From the repository root, run the following with Python 3 and a POSIX shell. The Lean check also requires the pinned Lean 4.30.0 toolchain and uses its bundled standard library.

```sh
(cd proof && shasum -a 256 -c SHA256SUMS && sh verify.sh)
(cd formalization && sh verify.sh)
```

Epoch asks for an expanded JSON array of all n+1 coefficients. The compact witness and mathematical proofs are available here. The expanded array has not been produced or submitted. Independent human peer review and Epoch acceptance are outstanding.
