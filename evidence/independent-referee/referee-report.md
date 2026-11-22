# Independent referee report: Dense Polynomials Whose Squares Are Sparse

**Review date:** 22 November 2026.

## Verdicts

**Mathematical proof: passes this independent audit.** The supplied manuscript and its exact compressed local-factor data give a specific finite polynomial in Z[x], prove that every coefficient from degree zero to its degree is nonzero, and establish the strict bound

\[
\operatorname{wt}(P^2)<\frac{501}{1000}\sqrt{\deg P}.
\]

I found no unsupported inference affecting the four mathematical issues requested. In particular, the endpoint corrections cover their joins with the interior, the residue restrictions control ordinary coefficients individually, and the finite evaluation arguments cover all relevant indices of all 128 local factors.

**Epoch submission compliance: fails the stated format requirement.** The supplied compressed masks and formulas are not the requested array of all coefficients. No expanded array was produced or submitted, and no run of Epoch's official verifier is established by this package. This conclusion is separate from the validity of the mathematical construction.

## Sources and audit scope

The sources are the uploaded `sparse-squares.pdf` and `epoch-sparse-squares-proof.zip`. The PDF inside the archive is byte-for-byte identical to the separately uploaded PDF. Every entry of the archive's SHA256SUMS manifest passed verification.

The official problem statement, accessed on 22 November 2026, explicitly requires a JSON array of plain integer literals in ascending degree order:

<https://epoch.ai/frontiermath/open-problems/polynomial-square-few-coefficients#prompt>

This report evaluates the construction on its own mathematical and computational evidence. It does not rely on the archive's prior-review statement, its saved verification log, or its finite-field fingerprints as a proof.

The audit included a fresh successful execution of the supplied `verify.sh`, source inspection of the verifiers, an independent derivation of the density and support arguments, independent regeneration of the 128 prime-selection records, and a separate checker that imports no proof-package code. That separate checker is included as `independent_check.py`; its captured output is `independent-check-results.json`.

The exact compressed-mask file audited has SHA-256:

```
f7e6aea9e98c666a7bcf7dcde5186d872234b047302d799918120331a54d4cee
```

## 1. A specific finite integer polynomial is defined

For each JSON record, let

\[
S_i(z)=\sum_j s_{i,j}z^j,
\qquad B_i(z)=\frac{S_i(z)}{(1-z)^4}.
\]

The records contain exact integer exponents and coefficients, with no duplicate exponents. For every record the four exact identities

\[
\sum_j s_{i,j}\binom{j}{r}=0\quad(r=0,1,2,3)
\]

hold. These are precisely the first four Taylor coefficients of S_i at 1, so (1-z)^4 divides S_i. Because this divisor is monic, the quotient is an integer polynomial. Its exact degree is max(supp S_i)-4, which is below 9v_i in every case.

The manuscript explicitly permits using these recorded polynomials directly in Sections 3-5. Thus their validity can be established without proving that a generator recreated the same numbers. This direct definition suffices for the claimed witness.

The manuscript also gives a deterministic local construction. Its linear system is invertible at y=1 because its matrix is vI-2J and v differs from twice the number of forbidden residues. The rank-two inverse is correct. Clearing denominators by the displayed multiplier gives integer sections of degree at most eight in y, hence local degree below 9v. The common multiplier is nonzero at 1 and has derivative zero there. The denominator and height estimates are conservative but valid.

The prime list is also unambiguous. The audit regenerated exactly the first 128 retained values under the prescribed scan of odd u and the prescribed test bases. The first prime is 461168601842738790401; the last is 116122253944001627422721, corresponding to u=6295.

The primality argument is deterministic. If q is a prime divisor of v=u2^64+1 and

\[
a^{u2^{63}}\equiv-1\pmod v,
\]

then the order of a modulo q divides u2^64 but not u2^63. Since u is odd, this order is divisible by 2^64, so q is at least 2^64+1. The verified inequality 2^64+1>sqrt(v) then rules out compositeness.

Every subsequent quantity is an explicitly defined integer, polynomial, finite product, or finite geometric sum. The ceiling defining b can be computed exactly as

\[
b=1+\operatorname{isqrt}((14t+3)N-1),\qquad t=128.
\]

There is no limiting process, unspecified choice, or need to expand the coefficient array in order to define P.

## 2. Local identities: the data really provide the needed hypotheses

Beyond divisibility, the separate checker verifies the local values and derivatives through sparse congruences, independently of the supplied verifier's residue interpolation.

Write beta=B_i(1). It is computed exactly as the fourth Taylor coefficient of S_i. The fifth Taylor coefficient is B_i'(1). Every record satisfies beta!=0 and B_i'(1)=0. The checker also verifies

\[
S_i+\beta(1-z)^4\equiv0\pmod{z^{v_i}-1},
\]

and

\[
(1-z)S_i'+4S_i\equiv0\pmod{z^{v_i}-1}.
\]

At each nontrivial v_i-th root of unity, these identities give B_i=-beta and B_i'=0. At 1 the already checked Taylor coefficients give B_i=beta and B_i'=0. Thus

\[
B_i^2\equiv\beta^2\pmod{(1-z^{v_i})^2}.
\]

The cyclic coefficient vector follows explicitly:

\[
[B_i]_{\mathrm{cyc}}(z)
=-\beta+\frac{2\beta}{v_i}I_{v_i}(z).
\]

Consequently the coefficient at residue 0 is beta(2/v_i-1), and at every other residue it is 2beta/v_i. They are nonzero. Integrality is consistent with the exact integer quotient and was checked directly for the common cyclic value.

The independent local height check uses

\[
\|B_i\|_1
\leq\sum_{j\leq d_i}|s_{i,j}|\binom{d_i-j+4}{4},
\qquad d_i=\deg B_i.
\]

This follows by expanding (1-z)^(-4), applying the triangle inequality to each coefficient, and summing the binomial coefficients. Every resulting upper bound is at most v_i^21.

The missing ordinary-square residues are certified separately, as described next. They do not follow merely from the cyclic-square congruence.

## 3. Why eight evaluations suffice, and why no indices are missed

Fix one local factor and write T=S^2=sum_j t_j z^j. The square is computed by ordinary sparse convolution. Then

\[
[z^k]B(z)^2
=\sum_{j\leq k}t_j\binom{k-j+7}{7}.
\]

Between successive exponents in supp(T), the set of active summands is constant. Each binomial expression is a polynomial in k of degree seven. Therefore on each such interval the coefficient is the restriction of one polynomial of degree at most seven.

Eight distinct zero values force that polynomial to vanish identically. If an interval has fewer than eight integer positions, checking every position suffices. This is an exact interpolation certificate, not a heuristic sampling argument.

The supplied `Coefficients` class correctly builds the polynomial

\[
7!\binom{k-j+7}{7}=\prod_{a=1}^{7}(k-j+a)
\]

and accumulates it as each numerator exponent becomes active. Its use of `bisect_right` includes the new numerator term at its own exponent.

The forbidden residues are exactly

\[
\ell=\lfloor v/3\rfloor+1,
\qquad h=\lfloor(2v-1)/3\rfloor.
\]

For each layer q from zero through floor(2d/v), where d=deg B, the verifier checks

\[
[qv+\ell,\ \min(qv+h+1,2d+1)).
\]

It splits this half-open interval at every interior exponent of S^2. Thus every forbidden exponent between 0 and 2d is covered. A knot at the left endpoint is active in that interval; a knot at the right endpoint is handled in the next one. Exponents above 2d are already zero because B^2 is a polynomial of degree 2d.

The separate checker uses a different acceptance test on long intervals: it checks that every coefficient of the degree-seven polynomial is exactly zero. It evaluates individual coefficients only on short intervals. This avoids relying solely on the same eight-point implementation.

### Exact coverage counts

| Quantity | Verified count |
|---|---:|
| Local factors | 128 |
| Nonzero terms in all S_i | 11,092 |
| Nonzero terms in all S_i^2 | 45,310 |
| Forbidden subintervals | 13,602 |
| Long subintervals, length at least eight | 2,044 |
| Short subintervals, tested in full | 11,558 |
| Supplied verifier's square-coefficient evaluations | 27,910 |
| Individual forbidden coefficients covered | 42,169,625,887,381,509,185,208,980 |
| Supplied verifier's cyclic/derivative evaluations | 3,324 |

For these specific data, every short forbidden subinterval has length one. Thus the separate checker verifies 2,044 polynomial identities plus 11,558 individual zero coefficients. The counts agree with the supplied verifier's 27,910 evaluations: 8*2044+11558=27910.

### Cyclic values, derivatives, and heights are also checked globally

For B=S/(1-z)^4, ordinary coefficients are piecewise cubic in their index. For a fixed layer q, substituting k=r+qv gives a piecewise cubic function of r; its possible formula changes occur at numerator exponents modulo v. The supplied verifier cuts at all these residues and also at 0,1,v. The cut at 1 isolates the exceptional residue zero.

The cyclic coefficient is therefore cubic on each interval, and the weighted sum sum_q(r+qv)B_{r+qv}, which checks the derivative condition, has degree at most four. Five distinct evaluations suffice; shorter intervals are tested in full. The layer range includes all coefficients, and the exact divisible-numerator formula gives zero beyond the degree. The supplied height verifier similarly uses exact Newton differences for a cubic on each interval. These are valid all-indices arguments.

## 4. Density of the auxiliary polynomial and endpoint corrections

The auxiliary V is dense. Put y=z^b. Its precursor K_b has sections

\[
(K_b)_0(y)=(b-1)^2-2y+y^2,
\]

and

\[
(K_b)_j(y)=-2(b+j)+2jy\quad(1\leq j<b).
\]

Multiplication by I_b(y) gives, for each nonzero residue j, an initial coefficient -2(b+j), b-1 interior coefficients -2b, and a final coefficient 2j. At residue zero the coefficients are (b-1)^2, (b-1)^2-2, b-2 entries equal to (b-1)^2-1, then -1 and 1. None vanishes for b>=3. Its degree is T=b(b+1).

The values and first derivatives of K_b at b-th roots of unity yield

\[
K_b^2=b^4+(1-z^b)^2Q_b,
\qquad\deg Q_b\leq2b.
\]

This proves wt(V^2)<=8b+2. Also (1-z)^2K_b has at most 3+6+9=18 terms, proving wt((1-z)^2V)<=18b.

Next let E=DH, F=EV(lambda x), and C=LA. The choice of lambda makes F dense and every cyclic coefficient of C nonzero. For 0<=j<=T, the coefficient of x^j in F, viewed as a polynomial in lambda, contains the nonzero term E_0 V_j lambda^j. For T<j<=T+deg E, it contains the nonzero term E_degE V_(j-deg E) lambda^(j-deg E). The condition T>=deg E puts the latter index in the range of V. Distinct V indices have distinct powers of lambda, so these particular terms cannot cancel symbolically.

The integer coefficients of these polynomials in lambda have magnitudes bounded by ||E||_1||V||_1. The specified lambda exceeds one plus this bound, so the Cauchy root bound rules out zero after specialization. Every cyclic coefficient of C has the form C_r(A)-2lambda C_(r-1)(A)+lambda^2 C_(r-2)(A), with a nonzero constant term and coefficients bounded by 2||A||_1. The same root bound applies.

Write d=deg C and w=deg F. The exact length conditions are

\[
d\leq w<MN,
\qquad \alpha>\|C\|_1.
\]

Using D I_M(x^N)=1-x^(MN), one has

\[
P=C I_M(x^N)+\alpha(1-x^{MN})F.
\]

Each coefficient of C I_M(x^N) is a sum of a subset of the coefficients of C, each appearing at most once. Its magnitude is therefore at most ||C||_1. Now split the full degree range into three exhaustive parts.

1. **0<=j<=w.** The correction coefficient is alpha F_j, a nonzero integer of magnitude at least alpha. It cannot be cancelled by the main summand.
2. **w<j<MN.** The correction is zero. Every exponent k of C congruent to j modulo N satisfies k<=d<=w<j and 0<=(j-k)/N<M. Consequently the coefficient is the full cyclic coefficient of C at that residue, already proved nonzero.
3. **MN<=j<=MN+w.** The correction is -alpha F_(j-MN), again too large to be cancelled by the main summand.

In particular j=w belongs to the first region, and j=MN belongs to the third. The adjacent interior indices are covered by the second. No join is omitted.

Finally deg(C I_M(x^N))=(M-1)N+d<MN+w. Thus the last correction coefficient survives, and

\[
\deg P=MN+w,
\qquad\operatorname{wt}(P)=MN+w+1.
\]

The audit checked both endpoint length inequalities using the actual degrees of the recorded local factors.

## 5. Why the residue count bounds ordinary coefficients

The local missing-residue statement is coefficientwise: B_i^2 has no ordinary term at any exponent whose residue belongs to R_i. In

\[
A^2=\prod_i B_i(x^{e_i})^2,
\]

all factors other than i have exponents divisible by v_i. Every exponent in the product must therefore lie, modulo v_i, in e_i((Z/v_i Z) minus R_i). Since e_i is invertible modulo v_i, the number of allowed residues modulo N is exactly

\[
s=\prod_i(v_i-|R_i|)\leq N(7/10)^t.
\]

This is a restriction on each ordinary exponent. It does not infer individual zeros from cancellation in a cyclic sum.

From the local double-root identities,

\[
A^2=c+(1-x^N)^2Q.
\]

Since zero is an allowed residue, subtracting c introduces no forbidden residue. In formal power series,

\[
Q=(A^2-c)\sum_{h\geq0}(h+1)x^{hN}.
\]

All shifts preserve the residue modulo N. Hence every coefficient in a forbidden residue remains zero individually. The quotient is already known to be a polynomial, so there is no convergence or infinite-support issue.

The degree bound deg Q<=18tN means each allowed residue contains at most 18t+1 candidate ordinary exponents. Therefore

\[
\operatorname{wt}(Q)\leq(18t+1)s
\leq(18t+1)N(7/10)^t.
\]

This is the required passage from residue restrictions to an ordinary support count.

## 6. The strict global inequality

The corrected square identity is exact:

\[
A_*^2=cL^2+D^2Q_*,
\]

\[
Q_*=L^2Q+2\alpha(AH)(LW)+\alpha^2D^2H^2W^2.
\]

The three support bounds give

\[
\operatorname{wt}(Q_*)\leq
5(18t+1)N(7/10)^t+18b10^{6t}+3(8b+2)9^t.
\]

Here wt(L^2)=5, wt(AH)<=10^(6t), wt(LW)<=18b, wt(D^2)=3, wt(H^2)<=9^t, and wt(W^2)<=8b+2. Nonzero scalar substitutions preserve support and products can only have fewer terms than the product of their support sizes.

For any certified integer upper bound q_* for wt(Q_*), define

\[
U=5(2M-1)+3q_*.
\]

The identity

\[
P^2=cL^2I_M(x^N)^2+(1-x^{MN})^2Q_*
\]

then implies wt(P^2)<=U.

Even the paper's coarser bounds establish the original threshold without delicate numerical evaluation: q_*<N/400 implies wt(P^2)<13N/400; M>=N/800 and (13/400)^2<1/800 imply wt(P^2)^2<MN<=deg P.

For the stronger constant, exact integer arithmetic independently confirms

\[
(1000U)^2<501^2 n,
\qquad n=\deg P.
\]

The integer n agrees with certificate.json and has 5,798 decimal digits. N has 2,901 decimal digits. The paper's conservative upper bound gives

\[
U/\sqrt n\approx0.5000000000000102882452.
\]

Using the actual sparse numerator weights, actual local degrees, and exact allowed-residue count gives a tighter upper bound U_data with

\[
0.500000000000000018\leq U_{\rm data}/\sqrt n
<0.500000000000000019.
\]

The separate checker verifies this last interval using integer square roots, and directly verifies (1000 U_data)^2<501^2 n. Both bounds prove the claimed 0.501 inequality.

Neither decimal is the actual support of the enormous square. In particular these bounds do not justify replacing 0.501 by the strict constant 0.5. Exact support enumeration is unnecessary for the claimed inequality, but it has not been performed.

## 7. Submission compliance and the requested gap assessment

Epoch's stated output is an expanded JSON array [a_0,a_1,...,a_n] containing plain, nonzero integer literals in ascending degree order. The local-mask JSON contains sparse numerator descriptions, not those coefficients. A generator, compressed definition, or correctness certificate is not that required array.

The degree has 5,798 digits, so the required array has more than 10^5797 entries. Its impractical size does not invalidate the finite mathematical definition, but the required file is absent. Thus the mathematical existence question is answered by this construction while the literal submission requirement remains unmet.

I found no first unsupported inference affecting the mathematical conclusions requested. Two weaker evidentiary claims should not be substituted for the actual proof: modular fingerprints do not prove equality of huge polynomials, and a few evaluations do not establish global zero coefficients without the degree and partition argument. This package supplies the requisite exact identities and partition argument, and the independent checks confirm the recorded data satisfy them.

## Reproduction

Extract the original proof ZIP and run, using Python 3.9 or newer:

```
python3 independent_check.py /path/to/epoch-sparse-squares
```

The separate checker reads the original JSON, period CSV, and degree certificate, but imports no Python module from that directory. It uses exact integer arithmetic and explicit failure checks, not Python assertions. It does not generate the global coefficient array or submit anything to Epoch.

The supplied program can separately be replayed from its original directory with `sh verify.sh`, with Python assertions enabled. The full original replay completed successfully in this audit.
