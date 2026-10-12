"""Synthetic CheckM2 DIAMOND DB whose KO labels follow real reference KO
profiles: for each test genome a CheckM2 reference genome with a similar CDS
count is chosen and its KO multiset is assigned to the genome's proteins."""
import os, sys, json, random, subprocess, pickle, warnings
import numpy as np, scipy.sparse as sp
warnings.filterwarnings('ignore')
random.seed(11)
out, cm2 = sys.argv[1], sys.argv[2]
fo = json.load(open(cm2 + '/data/feature_ordering.json'))
kos = fo['KO_Genes']
sc = pickle.load(open(cm2 + '/models/scaler.sav', 'rb'))
R = sp.load_npz(cm2 + '/data/min_ref_rsdata_v1.npz').tocsr()
def read_fasta(p):
    name, seq = None, []
    for l in open(p):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)
AA = 'ACDEFGHIKLMNPQRSTVWY'
refs, n = [], 0
for t in ['TEST1', 'TEST2', 'TEST3']:
    faa = '%s/%s.ref.faa' % (out, t)
    prots = [(nm, p.rstrip('*')) for nm, p in read_fasta(faa)]
    ncds = len(prots)
    # reference rows with a similar CDS count (column 21, unscaled)
    cds = (R[:, 21].toarray().ravel() - sc.min_[21]) / sc.scale_[21]
    cand = [i for i in range(R.shape[0]) if abs(cds[i] - ncds) / ncds < 0.15]
    row = random.choice(cand)
    r = R.getrow(row)
    pool = []
    for j, v in zip(r.indices, r.data):
        if 22 <= j < 22 + len(kos):
            c = int(round((v - sc.min_[j]) / sc.scale_[j]))
            pool += [kos[j - 22]] * c
    random.shuffle(pool)
    print(t, 'proteins', ncds, 'reference row', row, 'KO copies', len(pool))
    order = list(range(ncds)); random.shuffle(order)
    for pi, ko in zip(order, pool):
        p = prots[pi][1]
        if len(p) < 30: continue
        rate = random.choice([0.02, 0.1, 0.2, 0.35])
        m = ''.join(random.choice(AA) if random.random() < rate else ch for ch in p)
        if random.random() < 0.05:
            m = m[: int(len(m) * random.choice([0.4, 0.7]))]
        n += 1
        refs.append(('UniRef100_%07d~%s' % (n, ko), m))
for i in range(3000):
    n += 1
    refs.append(('UniRef100_%07d~%s' % (n, random.choice(kos)), ''.join(random.choice(AA) for _ in range(random.randint(80, 600)))))
os.makedirs(out + '/db/checkm2/CheckM2_database', exist_ok=True)
with open(out + '/ref.faa', 'w') as f:
    for nm, s in refs: f.write('>%s\n%s\n' % (nm, s))
subprocess.run(['diamond', 'makedb', '--in', out + '/ref.faa', '-d', out + '/db/checkm2/CheckM2_database/uniref100.KO.1', '--quiet'], check=True)
print('references', len(refs))
