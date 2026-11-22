"""Exact evaluation and interpolation check of the bounded-difference local masks."""
from argparse import ArgumentParser
from fractions import Fraction
from functools import reduce
from math import gcd, lcm


def convolve(left, right):
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                product[i + j] += a * b
    return product


def lagrange_basis(nodes):
    basis = []
    for i, node in enumerate(nodes):
        polynomial = [1]
        for j, other in enumerate(nodes):
            if j != i:
                polynomial = convolve(polynomial, [Fraction(-other, node-other), Fraction(1, node-other)])
        basis.append(polynomial)
    return basis


def local(p):
    middle = [j for j in range(p) if p < 3*j < 2*p]
    m = len(middle)
    free = set(range(p)).difference(middle)

    def rank_two(y):
        slope = Fraction(2, p)*(y - 1)
        w = [-2 - slope*j for j in middle]
        a = 1 + sum(w)/p
        b = sum(value*j for value, j in zip(w, middle))/p
        c = slope*m/p
        d = 1 + slope*sum(middle)/p
        return w, slope, a, b, c, d, a*d - b*c

    def determinant(y):
        return rank_two(Fraction(y))[-1]

    d1 = determinant(1)
    derivative = 2*d1 + (determinant(2) - determinant(0))/2

    def evaluate(point):
        y = Fraction(point)
        w, slope, a, b, c, d, det = rank_two(y)
        if not det:
            return None
        sections = {0: Fraction(p-2)}
        for j in free - {0}:
            if j <= p-m:
                sections[j] = -2*(1 - Fraction(j, p)*(y-1))
            else:
                sections[j] = Fraction(-2, y)*(1 - Fraction(j-p, p)*(y-1))
        rhs = []
        for r in middle:
            total = 0
            for i in free:
                j = (r-i) % p
                if j in free:
                    total += (y if i+j >= p else 1)*sections[i]*sections[j]
            rhs.append(-total/2)
        s1 = sum(v*h for v, h in zip(w, rhs))
        s2 = slope*sum(rhs)
        v1 = (d*s1 - b*s2)/det
        v2 = (a*s2 - c*s1)/det
        for j, h in zip(middle, rhs):
            sections[j] = h/p - (v1 + j*v2)/p**2
        clear = y**2*det*(d1 - derivative*(y-1))
        return [sections[j]*clear for j in range(p)]

    nodes, values = [], []
    point = 1
    while len(nodes) < 9:
        sample = evaluate(point)
        if sample is not None:
            nodes.append(point)
            values.append(sample)
        point += 1
    basis = lagrange_basis(nodes)
    coefficients = [sum(values[h][j]*basis[h][k] for h in range(9)) for k in range(9) for j in range(p)]
    scale = lcm(*(x.denominator for x in coefficients))
    integers = [int(x*scale) for x in coefficients]
    content = reduce(gcd, integers)
    primitive = [value//content for value in integers]
    while primitive and not primitive[-1]:
        primitive.pop()
    return primitive, middle


def checks(mask, p, middle):
    square = convolve(mask, mask)
    cyclic = [sum(mask[i::p]) for i in range(p)]
    derivative_sections = [sum(i*mask[i] for i in range(j, len(mask), p)) for j in range(p)]
    constant = sum(mask)**2
    defect = square.copy()
    defect[0] -= constant
    quotient = []
    for j, value in enumerate(defect):
        quotient.append(value + (2*quotient[j-p] if j >= p else 0) - (quotient[j-2*p] if j >= 2*p else 0))
    return {
        'degree': len(mask)-1 < 9*p,
        'cyclic': all(cyclic) and all(cyclic[s] == cyclic[1] for s in range(1, p)),
        'missing': all(not any(square[r::p]) for r in middle),
        'jet': not any(derivative_sections),
        'constant_square': bool(constant) and not any(quotient[-2*p:]),
        'fourth_difference': sum(map(bool, convolve(mask, [1,-4,6,-4,1]))) < 100000,
    }


def inspect(p, tamper=False):
    mask, middle = local(p)
    invariants = checks(mask, p, middle)
    assert all(invariants.values()), (p, invariants)
    fourth = sum(map(bool, convolve(mask, [1,-4,6,-4,1])))
    print(f'p={p}  m={len(middle)}  deg={len(mask)-1}  cyclic={p}/{p}  absent={middle[0]}..{middle[-1]}  jets=PASS  square_remainder=PASS  fourth={fourth}')
    if tamper:
        mask[1] += 1
        failed = [name for name, ok in checks(mask, p, middle).items() if not ok]
        assert {'missing', 'jet', 'constant_square'} <= set(failed)
        print(f'tamper p={p}: A[1] += 1; detected failures = {", ".join(failed)}')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('primes', type=int, nargs='*', default=[7,11,13,17,29,31,53,67,101,127,257])
    parser.add_argument('--tamper', action='store_true', help='change one coefficient of the first mask and verify rejection')
    args = parser.parse_args()
    for i, p in enumerate(args.primes):
        inspect(p, tamper=args.tamper and i == 0)
