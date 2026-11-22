"""Two-residue linear mask prototypes; polynomial coefficients are in ascending order."""
import math
from functools import reduce

def add(*items):
    vals=[sum(p[i] for p in items if i<len(p)) for i in range(max(map(len,items),default=0))]
    while vals and vals[-1]==0: vals.pop()
    return vals
def mul(a,b):
    if not a or not b: return []
    vals=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b): vals[i+j]+=u*v
    return add(vals)
def neg(a): return [-u for u in a]
def det(a,b,c,d): return add(mul(a,d),neg(mul(b,c)))
def template(p):
    R=(p//2,p//2+1)
    assert not set((u+v)%p for u in R for v in R).intersection(R)
    H=[[] if i in R else ([p*(p-2)] if not i else [-2*(p+i),2*i]) for i in range(p)]
    def shifted(f,k): return ([0]*k+f) if f else []
    def mixed(r,j):
        i=(r-j)%p
        return shifted([2*x for x in H[i]],int(i+j>=p))
    def free(r):
        return neg(add(*(shifted(mul(H[i],H[j]),(i+j)//p) for i in range(p) for j in range(p) if (i+j)%p==r)))
    (a,b),(c,d)=[[mixed(r,j) for j in R] for r in R]
    v,w=map(free,R)
    D=det(a,b,c,d)
    J=[sum(D)+sum(i*x for i,x in enumerate(D)),-sum(i*x for i,x in enumerate(D))]
    B=[mul(h,D) for h in H]
    B[R[0]]=det(v,b,w,d)
    B[R[1]]=det(a,v,c,w)
    B=[mul(h,J) for h in B]
    G=reduce(math.gcd,(v for h in B for v in h))
    B=[[v//G for v in h] for h in B]
    A=[0]*(p*max(map(len,B)))
    for i,h in enumerate(B): A[i::p]=h+[0]*(len(A[i::p])-len(h))
    A=add(A)
    square=mul(A,A); base=sum(A)
    for r in range(p):
        assert sum(square[r::p])==(base**2 if r==0 else 0)
        assert sum((r+p*k)*v for k,v in enumerate(square[r::p]))==0
        if r in R: assert not any(square[r::p])
    return B,A
if __name__=='__main__':
    for p in (5,7):
        B,A=template(p)
        print(p, B,'degree',len(A)-1,'missing', [i for i in range(p) if not any(mul(A,A)[i::p])])
