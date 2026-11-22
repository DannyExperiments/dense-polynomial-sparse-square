"""Independent exact verifier for the factors of one sparse-square witness.

The input is a list of sparse integer numerators. No local linear-system solver
is used. Between numerator knots, coefficients of T/(1-x)^d are polynomials
in the degree of degree at most d-1. Intersecting those intervals with the
forbidden sections allows an exact check independent of the size of a period.
"""
from bisect import bisect_right
import csv
import hashlib
import json
from math import comb, factorial, isqrt, prod
from pathlib import Path


class Coefficients:
    def __init__(self, terms, power):
        self.knots=sorted(terms)
        self.factor=factorial(power-1)
        self.polynomials=[]
        current=[0]*power
        for exponent in self.knots:
            basis=[1]
            for j in range(1,power):
                updated=[0]*(len(basis)+1)
                for k,coefficient in enumerate(basis):
                    updated[k]+=(j-exponent)*coefficient
                    updated[k+1]+=coefficient
                basis=updated
            current=[a+b*terms[exponent] for a,b in zip(current,basis)]
            self.polynomials.append(current)

    def __call__(self, exponent):
        i=bisect_right(self.knots,exponent)-1
        if i<0:return 0
        total=0
        for a in reversed(self.polynomials[i]):
            total=total*exponent+a
        value,rem=divmod(total,self.factor)
        assert not rem
        return value


def differences(values):
    while values:
        yield values[0]
        values=[b-a for a,b in zip(values,values[1:])]


def support_size(terms):
    return len(terms)


def convolution(terms):
    sq={}
    for a,x in terms.items():
        for b,y in terms.items():
            sq[a+b]=sq.get(a+b,0)+x*y
    return {k:v for k,v in sq.items() if v}


def test_local(p, terms):
    assert len(terms)<10**6 and all(isinstance(v,int) and v for v in terms.values())
    assert min(terms)>=0
    assert all(sum(v*comb(k,r) for k,v in terms.items())==0 for r in range(4))
    degree=max(terms)-4
    assert 0<=degree<9*p
    B=Coefficients(terms,4)

    # Bound the 1-norm from exact Newton interpolation on each run.
    breaks=sorted({0,degree+1,*[k for k in terms if k<=degree]})
    height=0
    for lo,hi in zip(breaks,breaks[1:]):
        size=hi-lo
        ds=list(differences([B(lo+j) for j in range(min(4,size))]))
        height+=size*sum(abs(v)*comb(size-1,j) for j,v in enumerate(ds))
    assert height<=p**21

    # Coefficients summed cyclically have degree at most 3 in the residue.
    # The derivative check has degree at most 4. Beyond deg(B) the formula
    # for B from the divisible numerator returns exactly zero.
    breaks=sorted({0,1,p,*[k%p for k in terms]})
    layers=1+degree//p
    samples=0
    common=sum(B(1+h*p) for h in range(layers))
    assert common
    for lo,hi in zip(breaks,breaks[1:]):
        for r in range(lo,min(lo+5,hi)):
            indices=[r+h*p for h in range(layers)]
            s=sum(B(i) for i in indices)
            assert (-2*s==(p-2)*common) if r==0 else s==common
            assert sum(i*B(i) for i in indices)==0
            samples+=1

    # Every ordinary coefficient of B² at a forbidden residue is zero.
    square=convolution(terms)
    B2=Coefficients(square,8)
    lo,hi=p//3+1,(2*p-1)//3
    checked=0
    for h in range(1+2*degree//p):
        start,stop=h*p+lo,min(h*p+hi+1,2*degree+1)
        if stop<=start:continue
        cuts=[start,*[k for k in B2.knots if start<k<stop],stop]
        for a,b in zip(cuts,cuts[1:]):
            for k in range(a,min(a+8,b)):
                assert B2(k)==0
                checked+=1
    return degree,len(terms),samples,checked


if __name__=='__main__':
    here=Path(__file__).parent
    with (here/'proth-periods.csv').open() as source:
        periods=[tuple(int(row[k]) for k in ('odd','prime','witness'))
                 for row in csv.DictReader(source)]
    primes=[p for _,p,_ in periods]
    assert len(primes)==len(set(primes))==128
    for odd,p,a in periods:
        assert odd%2 and 0<odd<2**64 and p==odd*2**64+1 and p>10**20
        assert pow(a,(p-1)//2,p)==p-1 and isqrt(p)<2**64+1

    raw=(here/'compact-local-masks.json').read_bytes()
    masks=json.loads(raw)
    assert [m['prime'] for m in masks]==primes
    results=[]
    for rec in masks:
        p=rec['prime']
        terms=dict(rec['terms'])
        assert len(terms)==len(rec['terms'])
        results.append(test_local(p,terms))
    N=prod(primes)
    t=len(primes)
    b=1+isqrt((14*t+3)*N-1)
    m=N//400
    degA=sum((N//p)*result[0] for p,result in zip(primes,results))
    degH=4*sum(N//p for p in primes)
    degE=N+degH
    degV=b*(b+1)
    degF=degE+degV
    assert max(degA+2,degE)<=degV and degF<m*N
    allowed=prod(p-((2*p-1)//3-p//3) for p in primes)
    q0=(1+(2*degA)//N)*allowed
    qnew=5*q0+18*b*prod(result[1] for result in results)+3*(8*b+2)*9**t
    assert 400*qnew<N
    ks=5*(2*m-1)+3*qnew
    n=m*N+degF
    assert ks**2<n
    print('SHA-256 of local certificates:',hashlib.sha256(raw).hexdigest())
    print('All',len(results),'prime periods certified')
    print('Range of actual sparse numerator weights:',min(r[1] for r in results),max(r[1] for r in results))
    print('All cyclic values and jets certified by polynomial identities from',sum(r[2] for r in results),'exact samples')
    print('All forbidden square coefficients certified by polynomial identities from',sum(r[3] for r in results),'exact samples')
    print('All 128 local heights and both endpoint length conditions certified')
    print('Strict global support bound certified:',ks**2<n)
    decimals=10**18
    ratio_floor=isqrt((ks*decimals)**2//n)
    print('Floor of 10^18 times the upper-bound ratio:',ratio_floor)
