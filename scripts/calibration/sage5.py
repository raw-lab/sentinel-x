import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, sys
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, 'sage2.py')).read().split("def discordant")[0])
def tau(L): return 0.06 + 0.025*np.sqrt(5000.0/L)
def discordant2(contigs, min_len=5000):
    cs=[c for c in contigs if len(c)>=min_len]
    if len(cs)<2: return 0.0
    P=np.array([tetra(c) for c in cs]); L=np.array([len(c) for c in cs],float)
    Pn=P/np.linalg.norm(P,axis=1,keepdims=True)
    D=1-Pn@Pn.T
    med=np.argmin(D@L)
    c=P[med].copy()
    assigned=None
    for _ in range(30):
        d=1-(P@c)/np.linalg.norm(P,axis=1)/np.linalg.norm(c)
        a=d<tau(L)
        if assigned is not None and (a==assigned).all(): break
        assigned=a
        c=(P[a]*L[a,None]).sum(0)/L[a].sum()
    return 1-L[assigned].sum()/L.sum()
for mean_kb in [8,15,40]:
    fps=[];errs=[]
    for f in G:
        for rep in range(3):
            fps.append(discordant2(fragment(list(G[f].values()), mean_kb=mean_kb)))
    for a,b in [('TEST1','TEST2'),('TEST1','TEST3'),('TEST3','TEST2'),('TEST3','TEST1'),('TEST2','TEST3'),('TEST2','TEST1')]:
        A=fragment(list(G[a].values()), mean_kb=mean_kb); B=fragment(list(G[b].values()), mean_kb=mean_kb)
        totA=sum(len(c) for c in A if len(c)>=5000)
        for frac in [0.05,0.15,0.3,0.45]:
            rng.shuffle(B); pick=[];s=0
            for cc in B:
                if len(cc)<5000: continue
                if s>=frac*totA: break
                pick.append(cc); s+=len(cc)
            d=discordant2(A+pick); truth=s/(totA+s); errs.append((d-truth,truth,a,b))
    e=np.array([x[0] for x in errs])
    print(f"mean_kb={mean_kb}: single-genome FP max={max(fps):.3f} mean={np.mean(fps):.4f} | mixture err mean={e.mean():+.4f} min={e.min():+.3f} max={e.max():+.3f}")
    for x in errs:
        if abs(x[0])>0.02: print("    ", x)
