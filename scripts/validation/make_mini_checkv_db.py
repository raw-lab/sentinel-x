"""Build a small CheckV database in the released layout from CheckV's own
test contigs, so real CheckV 1.1.1 and Sentinel-X can be run head to head.

  genome_db/checkv_reps.faa/.tsv/.dmnd, checkv_error.tsv
  hmm_db/checkv_hmms/{1,2,3}.hmm, checkv_hmms.tsv, genome_lengths.tsv

Families are built from pyrodigal-gv proteins of the test contigs:
substitution-only variant sets aligned without gaps and passed to hmmbuild.
Some contigs get microbial families on one or both ends (proviruses), some
references are near-identical (high-confidence AAI), some diverged.
"""
import os, sys, random, subprocess, gzip
import pyrodigal_gv

random.seed(20261008)
fna, out = sys.argv[1], sys.argv[2]
AA = 'ACDEFGHIKLMNPQRSTVWY'

def read_fasta(p):
    name, seq = None, []
    for l in open(p):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)

contigs = list(read_fasta(fna))
finder = pyrodigal_gv.ViralGeneFinder(meta=True)
genes = {}
for cid, seq in contigs:
    pg = finder.find_genes(seq.upper().encode())
    genes[cid] = [g.translate(include_stop=False) for g in pg]

def mutate(p, rate):
    return ''.join(random.choice(AA) if (c in AA and random.random() < rate) else c for c in p)

os.makedirs(out + '/genome_db', exist_ok=True)
os.makedirs(out + '/hmm_db/checkv_hmms', exist_ok=True)
tmp = out + '/build_tmp'; os.makedirs(tmp, exist_ok=True)

# ---- families
fam_rows, fam_files, nfam = [], {1: [], 2: [], 3: []}, 0
glen_rows = []
provirus_ends = {}
for ci, (cid, seq) in enumerate(contigs):
    ps = genes[cid]
    n = len(ps)
    if n == 0: continue
    # which genes become family seeds, and their category
    cats = ['none'] * n
    mode = ci % 7
    if mode == 0 and n >= 12:      # host on the left
        k = max(3, int(0.3 * n))
        for i in range(k): cats[i] = 'microbial'
        for i in range(k, n): cats[i] = 'viral' if random.random() < 0.5 else 'none'
    elif mode == 1 and n >= 15:    # host on both ends
        k = max(3, int(0.25 * n))
        for i in list(range(k)) + list(range(n - k, n)): cats[i] = 'microbial'
        for i in range(k, n - k): cats[i] = 'viral' if random.random() < 0.6 else 'none'
    elif mode == 2:                # no annotations at all
        pass
    else:
        for i in range(n): cats[i] = 'viral' if random.random() < 0.45 else 'none'
    for gi, cat in enumerate(cats):
        if cat == 'none' or len(ps[gi]) < 40: continue
        nfam += 1
        rate = random.choice([0.1, 0.2, 0.3, 0.4])
        members = [ps[gi]] + [mutate(ps[gi], rate) for _ in range(9)]
        name = ('VOG%05d' % nfam) if cat == 'viral' else ('PF%05d' % nfam)
        acc = None
        keyed_by_acc = (cat == 'microbial' and nfam % 3 == 0)
        if keyed_by_acc:
            acc = 'PF%05d.%d' % (nfam, random.randint(1, 30))
            name = 'Pfam_like_%d' % nfam
        sto = '%s/%s.sto' % (tmp, name)
        with open(sto, 'w') as f:
            f.write('# STOCKHOLM 1.0\n#=GF ID %s\n' % name)
            if acc: f.write('#=GF AC %s\n' % acc)
            for mi, m in enumerate(members):
                f.write('s%d %s\n' % (mi, m))
            f.write('//\n')
        hmm = '%s/%s.hmm' % (tmp, name)
        subprocess.run(['hmmbuild', '--amino', '-n', name, hmm, sto], check=True, stdout=subprocess.DEVNULL)
        fam_files[1 + nfam % 3].append(hmm)
        key = acc if keyed_by_acc else name
        cutoff = random.choice([10.0, 18.5, 25.0, 31.2, 40.0, 55.0])
        db = 'Pfam' if cat == 'microbial' else 'VOGDB'
        fam_rows.append((key, cat, cutoff, db))
        if cat == 'viral':
            ng = random.choice([4, 9, 10, 12, 25, 60])
            base = len(seq) * random.uniform(0.8, 1.6)
            cv = random.choice([0.0, 0.05, 0.15, 0.4])
            lens = sorted(max(500, int(random.gauss(base, base * cv))) for _ in range(ng))
            glen_rows.append((key, ng, cv, lens))
# a few HMMs that hit nothing and a tsv row without an HMM
for k in fam_files:
    with open('%s/hmm_db/checkv_hmms/%d.hmm' % (out, k), 'w') as f:
        for h in fam_files[k]: f.write(open(h).read())
with open(out + '/hmm_db/checkv_hmms.tsv', 'w') as f:
    f.write('hmm\tcategory\tscore_cutoff\tdatabase\n')
    for r in fam_rows: f.write('%s\t%s\t%s\t%s\n' % r)
    f.write('VOG99999\tviral\t20.0\tVOGDB\n')
with open(out + '/hmm_db/genome_lengths.tsv', 'w') as f:
    f.write('hmm\tgenomes\tcv\tlengths\n')
    for k, ng, cv, lens in glen_rows:
        f.write('%s\t%d\t%s\t%s\n' % (k, ng, cv, ','.join(map(str, lens))))

# ---- reference genomes (AAI)
reps, faa = [], []
nref = 0
for ci, (cid, seq) in enumerate(contigs):
    ps = genes[cid]
    if not ps or ci % 5 == 4: continue
    for rel in range(random.choice([1, 1, 2, 3])):
        nref += 1
        rid = ('DTR_%06d' if rel % 2 == 0 else 'GCA_%06d') % nref
        typ = 'circular' if rid.startswith('DTR') else 'genbank'
        rate = random.choice([0.01, 0.05, 0.15, 0.3, 0.45])
        frac = random.choice([1.0, 0.9, 0.7, 0.5])
        k = 0
        for p in ps:
            if len(p) < 30 or random.random() > frac: continue
            k += 1
            faa.append(('%s_%d' % (rid, k), mutate(p, rate)))
        if k == 0: continue
        length = int(len(seq) * random.choice([0.7, 0.95, 1.0, 1.1, 1.5]))
        reps.append((rid, typ, length))
with open(out + '/genome_db/checkv_reps.faa', 'w') as f:
    for n, s in faa: f.write('>%s\n%s\n' % (n, s))
with open(out + '/genome_db/checkv_reps.tsv', 'w') as f:
    f.write('checkv_id\thabitat\ttype\tlength\n')
    for rid, typ, length in reps: f.write('%s\tsynthetic\t%s\t%d\n' % (rid, typ, length))
subprocess.run(['diamond', 'makedb', '--in', out + '/genome_db/checkv_reps.faa', '-d', out + '/genome_db/checkv_reps', '--quiet'], check=True)

# ---- error table
with open(out + '/genome_db/checkv_error.tsv', 'w') as f:
    f.write('length\taai\tcov\tcount\terror\n')
    for L in [1000, 5000, 10000, 20000, 40000, 80000]:
        for aai in range(30, 101, 5):
            for cov in range(10, 101, 10):
                cnt = random.choice([0, 50, 101, 150, 500, 2000]) if aai < 60 else random.choice([150, 500, 2000])
                err = round(max(0.1, (100 - aai) * 0.4 + (100 - cov) * 0.15 + random.uniform(-1, 1)), 2)
                f.write('%d\t%d\t%d\t%d\t%s\n' % (L, aai, cov, cnt, err))
print('families', len(fam_rows), 'genome-length rows', len(glen_rows), 'references', len(reps), 'ref proteins', len(faa))
