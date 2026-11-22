"""Two-prime exact endpoint check, plus full modular square expansion.

Uses modular convolutions via carry-free integer packing. The full root
coefficients themselves are built and checked over the integers.
"""
from pathlib import Path
import sys
from math import isqrt
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'director'))
sys.path.insert(0,str(BASE/'sparse'))
from two_mask import template
from hermite import lift,periodize

def norm(a): return sum(abs(v) for v in a)
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: c[i+j]+=x*y
    return c
def add(a,b):
    c=a.copy()+[0]*max(0,len(b)-len(a))
    for i,y in enumerate(b): c[i]+=y
    while len(c)>1 and c[-1]==0: c.pop()
    return c
def scale(a,t): return [t*v for v in a]
def subst(a,d):
    b=[0]*(d*(len(a)-1)+1)
    b[::d]=a
    return b
def exactquot(f,step):
    # (1-z^N)^2 with nonzero constant; invert at the origin.
    q=[]
    for i,x in enumerate(f): q.append(x+(2*q[i-step] if i>=step else 0)-(q[i-2*step] if i>=2*step else 0))
    if any(q[-2*step:]): raise ValueError('nondivisible')
    return q[:-2*step]

def modmul(a,b,p):
    # Positive modular representatives; radices exceed maximal coefficient.
    wide=(min(len(a),len(b))*(p-1)**2).bit_length()//8+1
    def pack(vs): return int.from_bytes(b''.join((v%p).to_bytes(wide,'little') for v in vs),'little')
    out=pack(a)*pack(b)
    num=len(a)+len(b)-1
    raw=out.to_bytes(wide*num,'little')
    return [int.from_bytes(raw[i*wide:(i+1)*wide],'little')%p for i in range(num)]

def main():
    primes=(5,7)
    N=35
    U=H=[1]
    for p in primes:
        _,A=template(p)
        U=mul(U,subst(A,N//p))
        H=mul(H,subst([1,-4,6,-4,1],N//p))
    c=sum(U)**2
    D=subst([1,-1],N)
    Q=exactquot(add(mul(U,U),[-c]),N)
    KU=max(len(U)+1,len(mul(D,H))-1)
    k=max(3,isqrt(KU)+1)
    V=periodize(lift(k),k,k)
    assert all(V)
    T=len(V)-1
    assert T>=KU
    DH=mul(D,H)
    # Pick the first positive lambda surviving the exact polynomial tests.
    for lam in range(1,101):
        LV=[1,-2*lam,lam*lam]
        B=mul(U,LV)
        Vlam=[a*pow(lam,i) for i,a in enumerate(V)]
        E=mul(DH,Vlam)
        cyclic=[sum(B[i::N]) for i in range(N)]
        if all(cyclic) and all(E): break
    else: raise RuntimeError('small control needed a larger lambda')
    alpha=norm(B)+1
    W=len(E)-1
    M=W//N+5
    step=N*M
    A=add(B,scale(mul(D,E),alpha))
    P=periodize(A,N,M)
    I=subst([1]*M,N)
    P2=add(mul(I,B),scale(add(E,scale([0]*step+E,-1)),alpha))
    assert P==P2 and all(P) and len(P)==step+W+1
    wt=lambda u: sum(bool(x) for x in u)
    report={'primes':primes,'N':N,'degU':len(U)-1,'local_Q_weight':wt(Q),'wt_UH':wt(mul(U,H)),'k':k,'deg_V':T,'lambda':lam,'W':W,'M':M,'P_degree':len(P)-1,'nonzero_P_coefficients':wt(P),'max_P_coefficient_bits':max(abs(v).bit_length() for v in P)}
    # Direct exact square is practical at this small scale; six-digit λ suffices.
    q_actual=mul(P,P)
    R=add(add(mul(mul(LV,LV),Q),scale(mul(mul(U,H),mul(LV,Vlam)),2*alpha)),scale(mul(mul(mul(D,D),mul(H,H)),mul(Vlam,Vlam)),alpha*alpha))
    window=subst([1,-2,1],step)
    q_formula=add(scale(mul(mul(LV,LV),mul(I,I)),c),mul(window,R))
    assert q_actual==q_formula
    # Independent convolution algorithm on two primes.
    for prime in [65521,1000003]:
        assert modmul(P,P,prime)==[v%prime for v in q_formula]
    report.update({'R_weight':wt(R),'P_squared_weight':wt(q_actual),'strict_threshold':wt(q_actual)**2<len(P)-1,'exact_identity':True,'modular_crosschecks':[65521,1000003]})
    print(report)
if __name__=='__main__':main()
