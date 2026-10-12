import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, sys, csv
sys.path.insert(0, HERE)
from calib_lib import read_fasta, profiles
seqs=read_fasta(os.environ.get('CHECKV_TEST_FNA', 'test_sequences.fna'))
q={}
for row in csv.DictReader(open(os.environ.get('CHECKV_QUALITY', 'cv_out/quality_summary.tsv')), delimiter='\t'):
    q[row['contig_id']]=(row['provirus'], row['proviral_length'], row['contig_length'])
def best(Q):
    m=len(Q); cs=np.cumsum(Q,0); b=None
    for k in range(2,m-1):
        muL=cs[k-1]/k; muR=(cs[-1]-cs[k-1])/(m-k)
        bss=k*(m-k)/m*((muL-muR)**2).sum()
        if b is None or bss>b[0]:
            b=(bss,k,1-muL@muR/np.linalg.norm(muL)/np.linalg.norm(muR))
    return b
n=0;fl=0
for name,s in seqs.items():
    Q=profiles(s)
    if len(Q)<4: continue
    n+=1
    bss,k,d=best(Q)
    flag=d>=0.12; fl+=flag
    print(f"{name[:40]:40s} len={len(s):6d} m={len(Q):3d} D={d:.3f} k={k*5000:6d} {'FLAG' if flag else ''} provirus={q[name][0]} viral_len={q[name][1]}")
print(n,'assessed',fl,'flagged')
