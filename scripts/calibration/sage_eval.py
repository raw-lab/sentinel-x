import numpy as np
AA="ACDEFGHIKLMNPQRSTVWY"
def read_faa(p):
    out=[];buf=[]
    for line in open(p):
        if line.startswith('>'):
            if buf: out.append(''.join(buf).replace('*',''))
            buf=[]
        else: buf.append(line.strip())
    if buf: out.append(''.join(buf).replace('*',''))
    return out
def comp(p):
    v=np.array([p.count(a) for a in AA],float); s=v.sum(); return v/s if s else v
def current_index(prots):
    used=[p for p in prots if len(p)>=60]
    C=np.array([comp(p) for p in used]); mu=C.mean(0)
    cos=(C@mu)/np.linalg.norm(C,axis=1)/np.linalg.norm(mu); d=1-cos
    z=(d-d.mean())/d.std(); return 100*(z>2).mean()
def gmm1d(x, iters=200):
    # 2-component EM init by quantiles
    w=np.array([0.5,0.5]); mu=np.quantile(x,[0.25,0.75]); sd=np.array([x.std(),x.std()])/2
    for _ in range(iters):
        p=np.vstack([w[k]*np.exp(-0.5*((x-mu[k])/sd[k])**2)/sd[k] for k in range(2)])+1e-300
        r=p/p.sum(0)
        nk=r.sum(1); w=nk/len(x); mu=(r*x).sum(1)/nk; sd=np.sqrt((r*(x-mu[:,None])**2).sum(1)/nk)+1e-9
    ll2=np.log((np.vstack([w[k]*np.exp(-0.5*((x-mu[k])/sd[k])**2)/(sd[k]*np.sqrt(2*np.pi)) for k in range(2)])).sum(0)).sum()
    m1,s1=x.mean(),x.std(); ll1=(-0.5*((x-m1)/s1)**2-np.log(s1*np.sqrt(2*np.pi))).sum()
    n=len(x); dbic=(-2*ll1+2*np.log(n))-(-2*ll2+5*np.log(n))
    D=np.sqrt(2)*abs(mu[0]-mu[1])/np.sqrt(sd[0]**2+sd[1]**2)
    return w.min(), D, dbic
def garp(prots):
    used=[p for p in prots if len(p)>=100]
    C=np.array([comp(p) for p in used])
    G=C[:,[AA.index(a) for a in 'GARP']].sum(1); F=C[:,[AA.index(a) for a in 'FYMINK']].sum(1)
    return np.log((G+1e-3)/(F+1e-3))
def pc1(prots):
    used=[p for p in prots if len(p)>=100]
    C=np.array([comp(p) for p in used])+1e-3
    L=np.log(C); L=L-L.mean(1,keepdims=True)
    L=L-L.mean(0); U,S,Vt=np.linalg.svd(L,full_matrices=False); return L@Vt[0]
base=os.environ.get('TEST_PROTEINS', 'proteins').rstrip('/') + '/'   # TEST1/2/3.ref.faa
P={k:read_faa(base+k+'.ref.faa') for k in ['TEST1','TEST2','TEST3']}
rng=np.random.default_rng(1)
def mix(a,b,frac):
    A=P[a]; B=P[b]; nb=int(round(len(A)*frac/(1-frac)))
    nb=min(nb,len(B)); idx=rng.choice(len(B),nb,replace=False)
    return A+[B[i] for i in idx], nb/(len(A)+nb)
print("%-22s %6s %7s | %8s %6s %8s | %8s %6s %8s" % ("set","truth","SAGEidx","garp_w","D","dBIC","pc1_w","D","dBIC"))
for k in P:
    x=garp(P[k]); y=pc1(P[k])
    print("%-22s %6.2f %7.2f | %8.3f %6.2f %8.1f | %8.3f %6.2f %8.1f" % (k,0,current_index(P[k]),*gmm1d(x),*gmm1d(y)))
for a,b in [('TEST1','TEST2'),('TEST1','TEST3'),('TEST3','TEST2'),('TEST3','TEST1')]:
    for f in [0.05,0.1,0.2,0.35,0.5]:
        prots,truth=mix(a,b,f)
        x=garp(prots); y=pc1(prots)
        print("%-22s %6.2f %7.2f | %8.3f %6.2f %8.1f | %8.3f %6.2f %8.1f" % (f"{a}+{b}",truth,current_index(prots),*gmm1d(x),*gmm1d(y)))
