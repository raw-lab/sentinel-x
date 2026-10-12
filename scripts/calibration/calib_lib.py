import sys, itertools, random
import numpy as np
from pathlib import Path

def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        line=line.strip()
        if line.startswith('>'):
            if name: seqs[name]=''.join(buf).upper()
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line)
    if name: seqs[name]=''.join(buf).upper()
    return seqs

comp = str.maketrans('ACGT','TGCA')
kmers = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
canon = {}
idx = {}
for k in kmers:
    rc = k.translate(comp)[::-1]
    c = min(k, rc)
    if c not in idx: idx[c]=len(idx)
    canon[k]=idx[c]
D = len(idx)  # 136

def tetra(seq):
    v = np.zeros(D)
    for i in range(len(seq)-3):
        j = canon.get(seq[i:i+4])
        if j is not None: v[j]+=1
    s=v.sum()
    return v/s if s>0 else v

def profiles(seq, W=5000):
    m = len(seq)//W
    return np.array([tetra(seq[j*W:(j+1)*W]) for j in range(m)])

def two_seg(P, minw=2):
    m = len(P)
    mu = P.mean(0)
    tss = ((P-mu)**2).sum()
    best=(0,None,0)
    cs = np.cumsum(P,0)
    for k in range(minw, m-minw+1):
        muL = cs[k-1]/k; muR=(cs[-1]-cs[k-1])/(m-k)
        bss = k*(m-k)/m*((muL-muR)**2).sum()
        r = bss/tss if tss>0 else 0
        if r>best[0]:
            cosLR = muL@muR/np.linalg.norm(muL)/np.linalg.norm(muR)
            best=(r,k,1-cosLR)
    return best

def perm_p(P, obs, n=200, seed=1):
    rng=np.random.default_rng(seed)
    c=0
    for _ in range(n):
        Q=P[rng.permutation(len(P))]
        if two_seg(Q)[0]>=obs: c+=1
    return (c+1)/(n+1)

