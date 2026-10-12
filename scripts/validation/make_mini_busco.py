"""Mini BUSCO lineages (odb10, odb12, odb12.2 layouts) and a protein set
exercising single, duplicated, fragmented, very large, missing and
multi-BUSCO genes."""
import os, sys, random, subprocess, shutil
random.seed(3)
src, out = sys.argv[1], sys.argv[2]
AA = 'ACDEFGHIKLMNPQRSTVWY'
def read_fasta(p):
    name, seq = None, []
    for l in open(p):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)
def mutate(p, r): return ''.join(random.choice(AA) if random.random() < r else c for c in p)
prots = [(n, p.rstrip('*')) for n, p in read_fasta(src)]
cands = [x for x in prots if len(x[1]) >= 150]
random.shuffle(cands)
markers = cands[:60]
tmp = out + '/build'; os.makedirs(tmp, exist_ok=True)
ids = ['%dat2' % (1000 + i) for i in range(len(markers))]
for bid, (n, p) in zip(ids, markers):
    sto = f'{tmp}/{bid}.sto'
    with open(sto, 'w') as f:
        f.write('# STOCKHOLM 1.0\n')
        for k in range(10): f.write(f's{k} {p if k == 0 else mutate(p, 0.25)}\n')
        f.write('//\n')
    subprocess.run(['hmmbuild', '--amino', '-n', bid, f'{tmp}/{bid}.hmm', sto], check=True, stdout=subprocess.DEVNULL)
# self scores
with open(f'{tmp}/markers.faa', 'w') as f:
    for bid, (n, p) in zip(ids, markers): f.write(f'>{bid}_self\n{p}\n')
scores = {}
for bid in ids:
    r = subprocess.run(['hmmsearch', '--tblout', '/dev/stdout', '-o', '/dev/null', f'{tmp}/{bid}.hmm', f'{tmp}/markers.faa'], capture_output=True, text=True).stdout
    for l in r.splitlines():
        if not l.startswith('#') and l.split()[0] == bid + '_self':
            scores[bid] = float(l.split()[5])
# input proteins
inp = [(n, p) for n, p in prots]
byname = {n: i for i, (n, _) in enumerate(inp)}
plan = {}
for i, (bid, (n, p)) in enumerate(zip(ids, markers)):
    k = i % 12
    if k == 0: plan[bid] = 'missing'; inp[byname[n]] = (n, ''.join(random.choice(AA) for _ in range(len(p))))
    elif k == 1: plan[bid] = 'dup'; inp.append((n + '_dup', mutate(p, 0.08)))
    elif k == 2: plan[bid] = 'frag'; inp[byname[n]] = (n, p[: len(p) // 2])
    elif k == 3: plan[bid] = 'large'; inp[byname[n]] = (n, p + mutate(p, 0.9) + ''.join(random.choice(AA) for _ in range(len(p) // 2)))
    elif k == 4: plan[bid] = 'paralog'; inp.append((n + '_par', mutate(p, 0.35)))
    else: plan[bid] = 'single'
# a fusion protein matching two single-copy markers
a, b = markers[5][1], markers[6][1]
inp.append(('fusion_5_6', mutate(a, 0.05) + mutate(b, 0.05)))
os.makedirs(out + '/input', exist_ok=True)
with open(out + '/input/genome.faa', 'w') as f:
    for n, p in inp: f.write(f'>{n}\n{p}\n')
lens = {bid: len(p) for bid, (n, p) in zip(ids, markers)}
for name, kind in [('mini_odb10', 'odb10'), ('mini_odb12', 'odb12'), ('mini_odb12.2', 'odb12.2')]:
    d = f'{out}/lineages/{name}'
    os.makedirs(d + '/hmms', exist_ok=True)
    for bid in ids: shutil.copy(f'{tmp}/{bid}.hmm', f'{d}/hmms/{bid}.hmm')
    with open(d + '/scores_cutoff', 'w') as f:
        for bid in ids: f.write(f'{bid}\t{scores[bid] * 0.5:.2f}\n')
    if kind == 'odb10':
        with open(d + '/lengths_cutoff', 'w') as f:
            for j, bid in enumerate(ids): f.write(f'{bid}\t0\t{0 if j == 7 else lens[bid] * 0.05:.2f}\t{lens[bid]}\n')
    elif kind == 'odb12.2':
        with open(d + '/lengths_cutoff', 'w') as f:
            for bid in ids: f.write(f'{bid}\t{lens[bid]}\t{lens[bid] * 0.04:.2f}\n')
    with open(d + '/dataset.cfg', 'w') as f:
        f.write(f'name={name}\nspecies=test\ndomain=prokaryota\ncreation_date=2026-10-08\nnumber_of_BUSCOs={len(ids)}\nnumber_of_species=1\nprodigal_genetic_code=11\n')
print({k: list(plan.values()).count(k) for k in set(plan.values())}, 'input proteins', len(inp))
