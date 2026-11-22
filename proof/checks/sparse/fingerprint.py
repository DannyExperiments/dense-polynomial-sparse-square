"""Evaluate the exact compact polynomial at a finite-field point without expanding it."""
from fractions import Fraction as F
from math import comb, isqrt, prod
from csv import DictReader
from pathlib import Path
from local_sections import Evaluation, differences, clearing


def norm(v,m):
    v=F(v)
    return v.numerator*pow(v.denominator,-1,m)%m


def binomial_geometric(length, z, q):
    if z==1:
        return [comb(length,k+1)%q for k in range(4)]
    zN=pow(z,length,q)
    inv=pow(1-z,-1,q)
    sums=[(1-zN)*inv%q]
    for k in range(1,4):
        sums.append((z*sums[-1]-zN*comb(length,k))*inv%q)
    return sums


def local_value(p,z,y,q):
    e=Evaluation(p,y)
    jet=list(differences([clearing(p,F(k)) for k in range(1,6)]))
    derivative=sum((-1)**(j-1)*jet[j]/j for j in range(1,len(jet)))
    scalar=2*pow(p,11,q)*norm(jet[0]-derivative*(y-1),q)%q
    total=0
    bounds=e.piece_boundaries
    for lo,stop in zip(bounds,bounds[1:]):
        length=stop-lo
        vals=[norm(e.times_denominator(j),q) for j in range(lo,lo+min(4,length))]
        sums=binomial_geometric(length,z,q)
        total+=pow(z,lo,q)*sum(s*v for s,v in zip(sums,differences(vals)))
    return scalar*total%q


def geom(length,x,q):
    if x==1:return length%q
    return (1-pow(x,length,q))*pow(1-x,-1,q)%q


def evaluate(rows,z,q):
    if not 0<z<q:raise ValueError('choose a nonzero residue')
    t=len(rows)
    N=prod(rows)
    b=1+isqrt((14*t+3)*N-1)
    M=N//400
    y=pow(z,N,q)
    lam=(2+12*pow(16,t,q)*pow(b,4,q)+2*pow(N,21,q))%q
    alpha=(1+(1+lam)**2*pow(N,21,q))%q
    a=h=1
    for p in rows:
        s=pow(z,N//p,q)
        a=a*local_value(p,s,y,q)%q
        h=h*pow(1-s,4,q)%q
    w=lam*z%q
    wb=pow(w,b,q)
    K=(b*b - 2*pow(geom(b,w,q),2,q)*(b-(b-1)*w)+(1-wb)**2)%q
    V=K*geom(b,wb,q)%q
    anew=((1-w)**2*a + alpha*(1-y)**2*h*V)%q
    return anew*geom(M,y,q)%q


if __name__=='__main__':
    import sys
    with Path(__file__).with_name("proth-periods.csv").open() as handle:
        rows=[int(row["prime"]) for row in DictReader(handle)]
    q=1_000_000_007
    for point in map(int,sys.argv[1:] or [1,2,17]):
        print(q,point,evaluate(rows,point,q),flush=True)
