"""Compare CheckV output directories table by table (exact, then numeric)."""
import sys, math, os
a, b = sys.argv[1], sys.argv[2]
files = ['quality_summary.tsv', 'completeness.tsv', 'contamination.tsv', 'complete_genomes.tsv',
         'tmp/gene_features.tsv', 'tmp/aai.tsv', 'viruses.fna', 'proviruses.fna']
def load(p):
    if not os.path.exists(p): return None
    return open(p).read().splitlines()
tot_bad = 0
for f in files:
    A, B = load(os.path.join(a, f)), load(os.path.join(b, f))
    if A is None or B is None:
        print(f'{f:28s} missing in', 'A' if A is None else 'B'); continue
    if f.endswith('.fna'):
        same = A == B
        print(f'{f:28s} {"IDENTICAL" if same else "DIFFERENT"} ({len(A)//2} records vs {len(B)//2})')
        if not same: tot_bad += 1
        continue
    key = lambda l: tuple(l.split('\t')[:2]) if f.endswith('gene_features.tsv') or f.endswith('aai.tsv') else l.split('\t')[0]
    ha, hb = A[0], B[0]
    if ha != hb: print(f'{f}: header differs'); tot_bad += 1
    ra = {key(l): l for l in A[1:]}; rb = {key(l): l for l in B[1:]}
    exact = sum(1 for k in ra if k in rb and ra[k] == rb[k])
    diffs = []
    cols = ha.split('\t')
    for k in ra:
        if k not in rb: diffs.append((k, 'missing in B')); continue
        if ra[k] == rb[k]: continue
        fa, fb = ra[k].split('\t'), rb[k].split('\t')
        for c, (x, y) in enumerate(zip(fa, fb)):
            if x == y: continue
            try:
                if math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=1e-12):
                    diffs.append((k, f'{cols[c]}: {x} vs {y} (within 1e-9)')); continue
            except ValueError: pass
            diffs.append((k, f'{cols[c]}: {x} | {y}'))
    for k in rb:
        if k not in ra: diffs.append((k, 'missing in A'))
    print(f'{f:28s} rows A={len(ra)} B={len(rb)} identical={exact} differing={len(set(d[0] for d in diffs))}')
    for d in diffs[:12]: print('     ', d)
    if diffs: tot_bad += 1
print('tables with differences:', tot_bad)
