import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, sys
sys.path.insert(0, HERE)
from calib_lib import read_fasta, tetra
base=os.environ.get('TEST_BINS', 'bins').rstrip('/') + '/'   # TEST1/2/3.fna
G={}
for f in ['TEST1','TEST2','TEST3']:
    G[f]=read_fasta(base+f+'.fna')
rng=np.random.default_rng(3)
def fragment(seqs, mean_kb=40, min_len=5000):
    out=[]
    for s in seqs:
        i=0
        while i<len(s):
            L=int(np.clip(rng.lognormal(np.log(mean_kb*1000),0.8),2000,400000))
            out.append(s[i:i+L]); i+=L
    return out
def discordant(contigs, tau, min_len=5000):
    cs=[c for c in contigs if len(c)>=min_len]
    if len(cs)<2: return 0.0, len(cs)
    P=np.array([tetra(c) for c in cs]); L=np.array([len(c) for c in cs],float)
    Pn=P/np.linalg.norm(P,axis=1,keepdims=True)
    D=1-Pn@Pn.T
    n=len(cs); parent=list(range(n))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    for i in range(n):
        for j in range(i+1,n):
            if D[i,j]<tau: parent[find(i)]=find(j)
    comp={}
    for i in range(n): comp.setdefault(find(i),0.0); comp[find(i)]+=L[i]
    dom=max(comp.values())
    return 1-dom/L.sum(), n
res=[]
for tau in [0.04,0.06,0.08]:
  print("tau",tau)
  for f in G:
    for rep in range(3):
        cs=fragment(list(G[f].values()))
        d,n=discordant(cs,tau); print(f"  {f} frag rep{rep} contigs={n} discordant={d:.3f}")
  for a,b in [('TEST1','TEST2'),('TEST1','TEST3'),('TEST3','TEST2'),('TEST3','TEST1'),('TEST2','TEST3')]:
    A=fragment(list(G[a].values())); B=fragment(list(G[b].values()))
    totA=sum(map(len,A))
    for frac in [0.05,0.1,0.25]:
        rng.shuffle(B); pick=[];s=0
        for c in B:
            if s>=frac*totA: break
            pick.append(c); s+=len(c)
        d,n=discordant(A+pick,tau)
        truth=s/(totA+s)
        print(f"  {a}+{b} truth={truth:.3f} contigs={n} discordant={d:.3f}")
