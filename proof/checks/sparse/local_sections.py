"""Evaluate the nine integer coefficients of any local section, also for large primes."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import sys


def integer_sums(lower, upper):
    """Return sums of 1, i, i^2 for lower <= i <= upper."""
    if lower > upper:
        return 0, 0, 0
    count = upper-lower+1
    linear = (lower+upper)*count//2
    squares = (upper*(upper+1)*(2*upper+1) - (lower-1)*lower*(2*lower-1))//6
    return count, linear, squares


def differences(values):
    values = list(values)
    while values:
        yield values[0]
        values = [b-a for a,b in zip(values, values[1:])]


def free_intervals(prime, y):
    """Each tuple specifies an affine free section c + s*j on [a,b]."""
    low = prime//3+1
    high = (2*prime-1)//3
    m = high-low+1
    rate = 2*(y-1)/prime
    parts = [(0,0,F(prime-2),F(0)),
             (1, low-1, F(-2), rate),
             (high+1, prime-m, F(-2), rate),
             (prime-m+1, prime-1, (-2-prime*rate)/y, rate/y)]
    return [(a,b,c,s) for a,b,c,s in parts if a<=b]


def multiply_affine_rows(irow, jrow, degree):
    a,b,c,s = irow
    d,e,f,r = jrow
    lo,hi=max(a,degree-e),min(b,degree-d)
    h,k,l = integer_sums(lo,hi)
    return (c*(f+r*degree)*h + (s*(f+r*degree)-c*r)*k - s*r*l)


def matrix(prime, y):
    lower, upper = prime//3+1,(2*prime-1)//3
    count, linear, squares = integer_sums(lower, upper)
    rate=2*(y-1)/prime
    a=1+(-2*count-rate*linear)/prime
    b=(-2*linear-rate*squares)/prime
    c=rate*count/prime
    d=1+rate*linear/prime
    return a,b,c,d


def clearing(prime, y):
    a,b,c,d = matrix(prime,y)
    return y*y*(a*d-b*c)


class Evaluation:
    def __init__(self, prime, y):
        self.p,self.y=prime,F(y)
        self.lower,self.upper=prime//3+1,(2*prime-1)//3
        self.rows=free_intervals(prime,self.y)
        self.mat=matrix(prime,self.y)
        boundary={self.lower,self.upper+1}
        # All transitions of the convolution endpoints occur at these sums.
        for a,b,_,_ in self.rows:
            for d,e,_,_ in self.rows:
                for shift in (0,prime):
                    for ab in (a,b):
                        for de in (d,e):
                            for delta in (0,1):
                                r=ab+de-shift+delta
                                if self.lower<r<self.upper+1:
                                    boundary.add(r)
        bounds=sorted(boundary)
        self.piece_boundaries=sorted(set(bounds + [j for a,b,_,_ in self.rows for j in (a,b+1)]))
        total=weighted=F(0)
        for lo,stop in zip(bounds,bounds[1:]):
            length=stop-lo
            points=[self.rhs(lo+i) for i in range(min(5,length))]
            total+=sum(comb(length,j+1)*v for j,v in enumerate(differences(points)))
            weighted+=sum(comb(length,j+1)*v for j,v in enumerate(
                differences([(lo+i)*x for i,x in enumerate(points)])))
        rate=2*(self.y-1)/prime
        self.vg=(-2*total-rate*weighted,rate*total)

    @lru_cache(None)
    def rhs(self, r):
        return sum(multiply_affine_rows(a,b,r)+self.y*multiply_affine_rows(a,b,r+self.p)
                   for a in self.rows for b in self.rows)

    def times_denominator(self,j):
        p,y=self.p,self.y
        a,b,c,d=self.mat
        det=a*d-b*c
        if j<self.lower or j>self.upper:
            for lo,hi,c0,s in self.rows:
                if lo<=j<=hi:
                    return (c0+s*j)*y*y*det
            raise ValueError(j)
        g0,g1=self.vg
        projected=(d-j*c)*g0+(-b+j*a)*g1
        return -y*y*(det*self.rhs(j)/p-projected/(p*p))/2


@lru_cache(maxsize=8)
def context(prime):
    evals=[Evaluation(prime,k) for k in range(1,10)]
    jet=list(differences([clearing(prime,F(k)) for k in range(1,6)]))
    derivative=sum((-1)**(j-1)*jet[j]/j for j in range(1,len(jet)))
    return evals,jet[0],derivative


def section(prime, j):
    if not 0<=j<prime:
        raise ValueError("residue out of range")
    evals,d1,derivative=context(prime)
    scale=2*prime**11
    val=[scale*e.times_denominator(j)*(d1-derivative*(e.y-1)) for e in evals]
    polynomial=[F(0)]*9
    basis=[1]
    for i,delta in enumerate(differences(val)):
        for k,bi in enumerate(basis):
            polynomial[k]+=delta*bi/factorial(i)
        nxt=[0]*(len(basis)+1)
        for k,bi in enumerate(basis):
            nxt[k]-=(i+1)*bi
            nxt[k+1]+=bi
        basis=nxt
    assert all(a.denominator==1 for a in polynomial)
    return [int(a) for a in polynomial]


def test(prime):
    from hermite import mult
    coefs=[0]*(9*prime)
    for j in range(prime):
        b=section(prime,j)
        assert sum((j+prime*k)*c for k,c in enumerate(b))==0
        coefs[j::prime]=b
    sq=mult(coefs,coefs)
    lo,hi=prime//3+1,(2*prime-1)//3
    assert all(sq[i]==0 for i in range(len(sq)) if lo<=i%prime<=hi)
    cyc=[sum(coefs[i::prime]) for i in range(prime)]
    assert all(cyc) and len(set(cyc[1:]))==1 and cyc[0]*(-2)==cyc[1]*(prime-2)
    print(prime, 'cyclic square, omitted ordinary sections, all first jets verified')


if __name__=='__main__':
    if sys.argv[1]=='--test':
        for p in map(int,sys.argv[2:]):
            test(p)
    else:
        print(section(int(sys.argv[1]),int(sys.argv[2])))
