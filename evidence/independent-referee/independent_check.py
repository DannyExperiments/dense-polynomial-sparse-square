"""Independent audit of compressed local factors, using coefficient identities.

This file does not import any code from the supplied proof archive.
Long forbidden intervals are checked by zero polynomial coefficients, not by
evaluation at eight points. Short intervals are checked at every index.
"""
from argparse import ArgumentParser
from bisect import bisect_right
from collections import defaultdict
from csv import DictReader
from hashlib import sha256
from json import loads, dumps
from math import comb, isqrt, prod
from pathlib import Path
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(20000)
parser = ArgumentParser(description=__doc__)
parser.add_argument('proof_directory', nargs='?', type=Path,
                    default=Path(__file__).parent/'source/epoch-sparse-squares',
                    help='Directory extracted from epoch-sparse-squares-proof.zip')
ROOT = parser.parse_args().proof_directory.resolve()


def check(condition, label):
    if not condition:
        raise ArithmeticError(label)


def mul_linear(poly, constant):
    result = [0] * (len(poly) + 1)
    for i, a in enumerate(poly):
        result[i] += constant * a
        result[i + 1] += a
    return result


def cumulative_numerators(terms, order):
    """(order-1)! times each quotient-coefficient polynomial in the index."""
    current = [0] * order
    knots, polynomials = [], []
    for exponent, value in sorted(terms.items()):
        basis = [value]
        for shift in range(1, order):
            basis = mul_linear(basis, shift - exponent)
        current = [a + b for a, b in zip(current, basis)]
        knots.append(exponent)
        polynomials.append(current)
    return knots, polynomials


def at(poly, argument):
    total = 0
    for a in reversed(poly):
        total = total * argument + a
    return total


def local(record):
    p = record['prime']
    raw = record['terms']
    check(all(type(k) is int and type(a) is int and k >= 0 and a != 0
              for k, a in raw), 'integer sparse terms')
    terms = dict(raw)
    check(len(terms) == len(raw), 'unique exponents')
    check(len(terms) < 10**6, 'numerator support')
    degree = max(terms) - 4
    check(0 <= degree < 9*p, 'local degree')
    moments = [sum(a*comb(k, r) for k, a in raw) for r in range(6)]
    check(moments[:4] == [0]*4, 'fourfold exact divisibility')
    beta = moments[4]
    check(beta != 0 and moments[5] == 0, 'value and derivative at one')

    # S + B(1)*(1-z)^4 is zero modulo z^p-1. Hence B=-B(1)
    # at each nontrivial pth root, and its cyclic coefficients are explicit.
    residue = defaultdict(int)
    for k, a in raw:
        residue[k % p] += a
    for j in range(5):
        residue[j] += beta*(-1)**j*comb(4, j)
    check(not any(residue.values()), 'nontrivial root values')
    check((2*beta) % p == 0, 'integral common cyclic value')
    common = 2*beta//p
    check(common != 0 and common-beta != 0, 'every cyclic coefficient nonzero')

    # (1-z)S' + 4S = (1-z)^5 B', as a polynomial identity.
    derivative = defaultdict(int)
    for k, a in raw:
        derivative[(k-1) % p] += k*a
        derivative[k % p] += (4-k)*a
    check(not any(derivative.values()), 'derivatives at all nontrivial roots')

    # The coefficientwise triangle bound summed exactly by hockey-stick.
    height = sum(abs(a)*comb(degree-k+4, 4) for k, a in raw if k <= degree)
    check(height <= p**21, 'local height')

    # Independent square convolution, once per unordered pair.
    square = defaultdict(int)
    rows = sorted(terms.items())
    for i, (k, a) in enumerate(rows):
        square[2*k] += a*a
        for ell, b in rows[i+1:]:
            square[k+ell] += 2*a*b
    square = {k: a for k, a in square.items() if a}
    knots, polys = cumulative_numerators(square, 8)
    check(not any(polys[-1]), 'square quotient tail polynomial zero')
    low, high = p//3 + 1, (2*p-1)//3
    long_runs = short_runs = short_evaluations = covered = 0
    for layer in range(2*degree//p + 1):
        left = layer*p + low
        right = min(layer*p + high + 1, 2*degree + 1)
        if left >= right:
            continue
        cuts = [left] + [k for k in knots if left < k < right] + [right]
        for begin, end in zip(cuts, cuts[1:]):
            index = bisect_right(knots, begin) - 1
            polynomial = polys[index] if index >= 0 else [0]*8
            length = end - begin
            covered += length
            if length >= 8:
                check(not any(polynomial), 'identically zero long forbidden run')
                long_runs += 1
            else:
                for k in range(begin, end):
                    check(at(polynomial, k) == 0, 'short forbidden run')
                short_runs += 1
                short_evaluations += length
    return dict(prime=p, degree=degree, terms=len(terms), square_terms=len(square),
                long_runs=long_runs, short_runs=short_runs,
                short_evaluations=short_evaluations, covered=covered)


def main():
    raw = (ROOT/'checks/sparse/compact-local-masks.json').read_bytes()
    records = loads(raw)
    check(len(records) == 128, '128 factors')
    with (ROOT/'proth-periods.csv').open() as file:
        rows = [{k: int(v) for k, v in row.items()} for row in DictReader(file)]
    primes = [row['prime'] for row in rows]
    check([record['prime'] for record in records] == primes, 'matching prime order')
    check(len(set(primes)) == 128, 'distinct primes')
    for row in rows:
        u, p, witness = row['odd'], row['prime'], row['witness']
        check(u % 2 == 1 and 0 < u < 2**64, 'Proth odd part')
        check(p == u*2**64+1 and p > 10**20, 'Proth form and size')
        check(pow(witness, (p-1)//2, p) == p-1, 'Proth witness')
        check((2**64+1)**2 > p, 'Proth prime divisor bound')
    results = [local(record) for record in records]
    N = prod(primes)
    t = len(primes)
    b = isqrt((14*t+3)*N-1)+1
    M = N//400
    dA = sum(r['degree']*(N//r['prime']) for r in results)
    dH = 4*sum(N//p for p in primes)
    dE = N+dH
    T = b*(b+1)
    w = T+dE
    check(T >= max(dA+2, dE), 'patch covers start of interior')
    check(w < M*N, 'disjoint endpoint patches')
    n = M*N+w
    allowed = prod(p-((2*p-1)//3-p//3) for p in primes)
    q0 = (1+2*dA//N)*allowed
    qstar = 5*q0+18*b*prod(r['terms'] for r in results)+3*(8*b+2)*9**t
    U = 5*(2*M-1)+3*qstar
    check(400*qstar < N, 'repaired quotient support')
    check(U*U < n, 'strict Epoch inequality')
    check((1000*U)**2 < 501**2*n, 'strict stronger inequality')
    certificate = loads((ROOT/'certificate.json').read_text())
    check(int(certificate['degree']) == n, 'recorded degree equality')
    scale = 10**18
    ratio_floor = isqrt((scale*U)**2//n)
    summary = dict(
        status='all independent checks passed',
        mask_sha256=sha256(raw).hexdigest(),
        factors=t,
        numerator_terms=sum(r['terms'] for r in results),
        squared_numerator_terms=sum(r['square_terms'] for r in results),
        long_forbidden_runs_tested_as_polynomial_identities=sum(r['long_runs'] for r in results),
        short_forbidden_runs_tested_exhaustively=sum(r['short_runs'] for r in results),
        short_run_coefficient_evaluations=sum(r['short_evaluations'] for r in results),
        forbidden_ordinary_coefficients_covered=sum(r['covered'] for r in results),
        degree_digits=len(str(n)),
        period_digits=len(str(N)),
        degree_prefix=str(n)[:32],
        degree_suffix=str(n)[-32:],
        floor_of_10_to_18_times_U_over_sqrt_n=ratio_floor,
        strict_501_over_1000_bound=True,
    )
    print(dumps(summary, indent=2))


if __name__ == '__main__':
    main()
