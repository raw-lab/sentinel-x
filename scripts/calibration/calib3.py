import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, sys
sys.path.insert(0, HERE)
from calib_lib import read_fasta, profiles
base=os.environ.get('TEST_BINS', 'bins').rstrip('/') + '/'   # TEST1/2/3.fna
P={}
for f in ['TEST1','TEST2','TEST3']:
    s=max(read_fasta(base+f+'.fna').values(), key=len)
    P[f]=profiles(s)

def stats(Q, minw=2):
    m=len(Q); mu=Q.mean(0); tss=((Q-mu)**2).sum()
    cs=np.cumsum(Q,0); best=None
    for k in range(minw, m-minw+1):
        muL=cs[k-1]/k; muR=(cs[-1]-cs[k-1])/(m-k)
        bss=k*(m-k)/m*((muL-muR)**2).sum()
        if best is None or bss>best[0]:
            d=1-muL@muR/np.linalg.norm(muL)/np.linalg.norm(muR)
            best=(bss,k,d)
    bss,k,d=best
    r2=bss/tss; wss=tss-bss
    F=bss/(wss/(m-2)) if wss>0 else 1e9
    return d,r2,F,k

def old_rule(Q):
    m=len(Q)
    mu=Q.mean(0)
    cos=lambda a,b: a@b/np.linalg.norm(a)/np.linalg.norm(b)
    dists=np.array([1-cos(w,mu) for w in Q]); mean=dists.mean(); sd=dists.std()
    z=(dists-mean)/sd if sd>1e-12 else np.zeros(m)
    maxz=max(0,z.max()); of=(z>2).mean()
    shifts=np.array([1-cos(Q[j-1],Q[j]) for j in range(1,m)])
    best=shifts.max(); bg=(shifts.sum()-best)/(len(shifts)-1) if len(shifts)>1 else 0
    seam = m>=3 and best>0.08 and best>bg*4+0.02
    broad = m>=8 and (of>0.15 or maxz>=4.0)
    return seam or broad

ctrl=[]  # (L, d, r2, F, old)
for L in [4,5,6,8,10,15,20,30,40,60,100,200]:
    for f,Q in P.items():
        step=1 if L<=10 else max(1,L//4)
        for o in range(0,len(Q)-L+1,step):
            R=Q[o:o+L]; d,r2,F,k=stats(R); ctrl.append((L,d,r2,F,old_rule(R)))
whole=[(f,stats(Q),old_rule(Q)) for f,Q in P.items()]
chim=[]
rng=np.random.default_rng(7); names=list(P)
for (a,b) in [(2,2),(2,3),(2,4),(2,8),(3,7),(5,5),(2,18),(4,16),(4,36),(10,30),(10,90),(2,98),(4,196),(20,180)]:
    for x in names:
        for y in names:
            if x==y: continue
            for _ in range(15):
                if len(P[x])<=a or len(P[y])<=b: continue
                oa=rng.integers(0,len(P[x])-a); ob=rng.integers(0,len(P[y])-b)
                Q=np.vstack([P[x][oa:oa+a],P[y][ob:ob+b]])
                d,r2,F,k=stats(Q); chim.append(((a,b),d,r2,F,old_rule(Q)))
np.save('ctrl.npy', np.array([(c[0],c[1],c[2],c[3],c[4]) for c in ctrl]))
np.save('chim.npy', np.array([(c[0][0],c[0][1],c[1],c[2],c[3],c[4]) for c in chim]))
print("whole chromosomes:", [(f, round(s[0],4), round(s[1],3), round(s[2],1), o) for f,s,o in whole])
C=np.array([(c[0],c[1],c[2],c[3],c[4]) for c in ctrl]); H=np.array([(c[0][0]+c[0][1],min(c[0]),c[1],c[2],c[3],c[4]) for c in chim])
print("old rule: control FP %.3f  chimera TPR %.3f" % (C[:,4].mean(), H[:,5].mean()))
for dt in [0.06,0.08,0.10,0.12]:
    for rt in [0.0,0.2,0.3,0.4,0.5]:
        fp=((C[:,1]>=dt)&(C[:,2]>=rt)).mean(); tp=((H[:,2]>=dt)&(H[:,3]>=rt)).mean()
        print(f"D>={dt:.2f} R2>={rt:.1f}: FP {fp:.4f} TPR {tp:.3f}")
