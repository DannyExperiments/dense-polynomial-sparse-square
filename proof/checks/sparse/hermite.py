"""Exact cyclic square roots from Hermite interpolation. No third-party dependencies."""
def mult(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if x:
            for j, y in enumerate(right):
                result[i + j] += x * y
    return result


def lift(period):
    n = period
    a = [0] * (2*n + 1)
    for j in range(n):
        a[j] = -2*(n+j)
        a[n+j] = 2*j
    a[0] += n*n + 1
    a[n] -= 2
    a[2*n] += 1
    return a


def periodize(a, n, m):
    p = [0] * (len(a) + (m-1)*n)
    for j in range(m):
        for i, value in enumerate(a):
            p[n*j+i] += value
    return p


def test(n, m):
    a = lift(n)
    p = periodize(a, n, m)
    q = mult(p, p)
    assert all(p)
    support = sum(x!=0 for x in q)
    assert support <= 2*m-1 + 3*(2*n+1)
    print('N M degree support bound', n,m,len(p)-1,support,2*m-1+3*(2*n+1))


if __name__=='__main__':
    for n in (3,4,5,7,13,17,31):
        test(n,3*n)
