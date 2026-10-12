import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, sys
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, 'sage5.py')).read().split("for mean_kb in")[0])
rng2=np.random.default_rng(11)
out=[]
for mean_kb,(a,b),frac in [(15,('TEST1','TEST3'),0.15),(8,('TEST3','TEST2'),0.3),(40,('TEST2','TEST1'),0.1)]:
    A=fragment(list(G[a].values()), mean_kb=mean_kb); B=fragment(list(G[b].values()), mean_kb=mean_kb)
    totA=sum(len(c) for c in A if len(c)>=5000); pick=[]; s=0
    for c in B:
        if len(c)<5000: continue
        if s>=frac*totA: break
        pick.append(c); s+=len(c)
    contigs=A+pick
    name=f"sim_{a}_{b}_{mean_kb}"
    with open(f'sim_bins/{name}.fna','w') as f:
        for i,c in enumerate(contigs): f.write(f">{name}_c{i}\n{c}\n")
    print(name, "python discordant %.4f" % (100*discordant2(contigs)), "truth %.4f" % (100*s/(totA+s)))
