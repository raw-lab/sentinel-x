"""A mini Pfam-like HMM library (NC cutoffs, names from viralVerify's
classifier table) from the mini CheckV families."""
import sys, os, random, glob, math
random.seed(5)
famdir, table, cvtsv, out = sys.argv[1:5]
rows = [l.split() for l in open(table)]
viral = sorted(rows, key=lambda r: math.log(float(r[6])) - math.log(float(r[7])), reverse=True)[:2500]
chrom = sorted(rows, key=lambda r: math.log(float(r[5])) - math.log(float(r[4])), reverse=True)[:2500]
plas = sorted(rows, key=lambda r: math.log(float(r[4])) - math.log(float(r[5])), reverse=True)[:2500]
cat = {}
for l in open(cvtsv).read().splitlines()[1:]:
    f = l.split('\t'); cat[f[0]] = f[1]
used = set()
def pick(pool):
    for _ in range(100000):
        r = random.choice(pool)
        if r[0] not in used:
            used.add(r[0]); return r[0]
    raise SystemExit('name pool exhausted')
with open(out, 'w') as o:
    for hmm in sorted(glob.glob(famdir + '/*.hmm')):
        txt = open(hmm).read()
        name = txt.split('NAME  ', 1)[1].split('\n', 1)[0].strip()
        acc = txt.split('ACC   ', 1)[1].split('\n', 1)[0].strip() if '\nACC   ' in txt else None
        c = cat.get(acc or name, cat.get(name, 'viral'))
        new = pick(viral if c == 'viral' else (plas if random.random() < 0.3 else chrom))
        txt = txt.replace('NAME  ' + name + '\n', 'NAME  ' + new + '\n', 1)
        nc = random.choice([15.0, 20.0, 25.0])
        txt = txt.replace('\nSTATS LOCAL MSV', f'\nNC    {nc:.2f} {nc - 3:.2f};\nSTATS LOCAL MSV', 1)
        o.write(txt)
print('HMMs', len(used))
