#!/usr/bin/env python3
"""Supplementary exact checks for the revised Glazer extractor.

Run from any directory with Python 3.10 or newer and no -O flag.
Uses only the standard library and the adjacent verifier. These checks
supplement, rather than replace, the uniform proofs in the manuscript.
"""
from pathlib import Path
import sys,json
from fractions import Fraction
from math import comb
from itertools import product
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
import verify_glazer_critical as v
if not __debug__:
    raise RuntimeError('Run without -O: assertions are part of this checker')
cert=json.loads((root/'certificates/finite_blocks.json').read_text())
for block in cert['blocks']:
    v.certificate_check(block)
print('Independent finite-witness checker: all six blocks passed')
for block in cert['blocks']:
    N=block['N']
    if N>16: continue
    a,R,f=block['a'],block['R'],block['f']
    signed=[0]*(2*N+1)
    count=0
    for i in range(N):
        Ri=R[i]
        ranks=[0]*(Ri+1)
        signs={}
        for bits in product('01',repeat=Ri):
            prefix=''.join(bits); r=prefix.count('0')
            e=(-1)**(i+r)*f[i][r]
            plus=(comb(Ri,r)+e)//2
            signs[prefix]=1 if ranks[r]<plus else -1
            ranks[r]+=1
        for bits in product('01',repeat=N-1-i):
            tail=''.join(bits)
            sign=signs[tail[:Ri]]
            z=N+i+tail.count('0')
            signed[z]+=sign
            signed[2*N-z]-=sign
            count+=2
    assert count==2*(2**N-1) and all(c==0 for c in signed)
    print(f'Deterministic prefix realization N={N}: {count} words, all signed type counts zero')
count=0
for n in range(33):
    for eta in [Fraction(1,2),Fraction(3,5),Fraction(7,10),Fraction(9,10)]:
        rho=2*(1-eta)
        target=[comb(n,k)*eta**k*(1-eta)**(n-k) for k in range(n+1)]
        mixture=[Fraction(0)]*(n+1)
        for m in range(n+1):
            weight=comb(n,m)*rho**m*(1-rho)**(n-m)
            for k in range(m+1):
                mixture[n-m+k]+=weight*Fraction(comb(m,k),2**m)
        assert target==mixture
        assert max(target)**2*(n+1)*rho<=1
        count+=1
print(f'Fair-coin mixture identity and modal square bound: {count} exact rational tests')
for N in range(1,65):
    assert all(comb(N,k)%2==0 for k in range(1,N)) == (N&(N-1)==0)
print('Dyadic parity criterion: N=1,...,64 checked')
# The large-block polynomial-identity path must reject a nonzero polynomial.
a,R=v.profile(256)
rows=[[0]*(ri+1) for ri in R]
rows[-1][0]=1
try:
    v.polynomial_identity(256,a,rows,1)
except AssertionError:
    print('Large-block radix checker rejects a deliberate nonzero polynomial')
else:
    raise AssertionError('checker accepted a false identity')
print('All independent supplementary checks passed.')
