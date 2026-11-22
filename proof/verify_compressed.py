"""Check sparse numerators and evaluate the witness without calling the local solver."""
from collections import defaultdict
from csv import DictReader
from json import loads
from math import comb, isqrt, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIELD = 1_000_000_007


def moment(terms, order):
    return sum(value * comb(index, order) for index, value in terms)


def validate_mask(prime, terms):
    positions = [k for k, _ in terms]
    assert positions == sorted(set(positions)) and positions[0] >= 0
    assert all(isinstance(v, int) and v for _, v in terms)
    assert 0 < len(terms) <= 10**6
    degree = positions[-1] - 4
    assert degree < 9*prime
    assert all(moment(terms, j) == 0 for j in range(4))
    at_one = moment(terms, 4)
    assert at_one and moment(terms, 5) == 0
    # The numerator T gives B=T/(1-z)^4, an integer polynomial.
    # B(1)=T^(4)(1)/4!. The next Taylor coefficient is zero.
    # At all other pth roots, B has value -B(1) and derivative zero.
    residues = defaultdict(int)
    differentiations = defaultdict(int)
    for k, value in terms:
        residues[k % prime] += value
        differentiations[(k-1) % prime] += k*value
        differentiations[k % prime] += (4-k)*value
    for j in range(5):
        residues[j] += at_one * (-1)**j * comb(4, j)
    assert not any(residues.values()) and not any(differentiations.values())
    # Coefficientwise bound from expanding (1-z)^-4, then summing up to deg B.
    norm_bound = sum(abs(value)*comb(degree-k+4, 4)
                     for k, value in terms if k <= degree)
    assert norm_bound <= prime**21
    return at_one


def geometric(x, terms):
    if x == 1:
        return terms % FIELD
    return (pow(x, terms, FIELD)-1) * pow(x-1, -1, FIELD) % FIELD


def evaluate(point, rows, period, b, lam, alpha):
    main = correction = 1
    for p, terms, at_one in rows:
        z = pow(point, period//p, FIELD)
        if z == 1:
            value = at_one
        else:
            numerator = sum(v*pow(z, k, FIELD) for k, v in terms)
            value = numerator*pow(pow(1-z, 4, FIELD), -1, FIELD)
        main = main*value % FIELD
        correction = correction*pow(1-z, 4, FIELD) % FIELD
    w = point*lam % FIELD
    wb = pow(w, b, FIELD)
    aux = (b*b - 2*geometric(w, b)**2*(b-(b-1)*w) + (1-wb)**2)
    aux = aux*geometric(wb, b) % FIELD
    y = pow(point, period, FIELD)
    lifted = (1-w)**2*main + alpha*(1-y)**2*correction*aux
    return lifted*geometric(y, period//400) % FIELD


def main():
    masks = loads((HERE/'checks/sparse/compact-local-masks.json').read_text())
    with (HERE/'proth-periods.csv').open() as stream:
        primes = [int(row['prime']) for row in DictReader(stream)]
    assert primes == [mask['prime'] for mask in masks]
    rows = [(p, mask['terms'], validate_mask(p, mask['terms']))
            for p, mask in zip(primes, masks)]
    period = prod(primes)
    b = isqrt((14*len(primes)+3)*period-1) + 1
    a = pow(period, 21, FIELD)
    lam = (2 + 12*pow(16, len(primes), FIELD)*pow(b, 4, FIELD) + 2*a) % FIELD
    alpha = (1+(1+lam)**2*a) % FIELD
    for line in (HERE/'fingerprints.txt').read_text().splitlines():
        field, point, expected = map(int, line.split())
        assert field == FIELD
        assert evaluate(point, rows, period, b, lam, alpha) == expected
    print("All 128 sparse numerators pass degree, integrality, value, derivative and height checks.")
    print("All four full-witness finite-field evaluations agree with the separate section evaluator.")


if __name__ == '__main__':
    main()
