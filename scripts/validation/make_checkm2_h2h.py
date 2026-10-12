"""Bins + a synthetic KO-labelled DIAMOND database (CheckM2 header format
`UniRef100_X~K00001`) for a CheckM2 vs Sentinel head-to-head."""
import os, sys, json, random, subprocess
random.seed(7)
src, out, fo = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out + '/bins', exist_ok=True)
def read_fasta(p):
    name, seq = None, []
    for l in open(p):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)
def write(p, recs):
    with open(p, 'w') as f:
        for n, s in recs:
            f.write('>%s\n' % n)
            for i in range(0, len(s), 80): f.write(s[i:i+80] + '\n')
g = {t: list(read_fasta('%s/%s.tst' % (src, t))) for t in ['TEST1', 'TEST2', 'TEST3']}
write(out + '/bins/TEST1.fna', g['TEST1'])
write(out + '/bins/TEST2.fna', g['TEST2'])
write(out + '/bins/TEST3.fna', g['TEST3'])
# fragmented partial bin: 50 kb pieces of TEST1's chromosome, ~40 % kept
c = g['TEST1'][0][1]
pieces = [(f'frag{i}', c[i:i+50000]) for i in range(0, len(c), 50000) if random.random() < 0.4]
write(out + '/bins/partial.fna', pieces)
# contaminated bin: TEST2 + TEST3
write(out + '/bins/mixed.fna', [('a_' + n, s) for n, s in g['TEST2']] + [('b_' + n, s) for n, s in g['TEST3']])
# small bin (< 100 kb -> Prodigal meta mode)
write(out + '/bins/small.fna', [('s1', g['TEST3'][0][1][200000:280000])])

# proteins for the reference: prodigal single on the full genomes
kos = json.load(open(fo))['KO_Genes']
AA = 'ACDEFGHIKLMNPQRSTVWY'
refs = []
n = 0
for t in ['TEST1', 'TEST2', 'TEST3']:
    faa = out + '/%s.ref.faa' % t
    subprocess.run(['prodigal', '-p', 'single', '-q', '-a', faa, '-i', out + '/bins/%s.fna' % t, '-o', '/dev/null'], check=True)
    for name, p in read_fasta(faa):
        p = p.rstrip('*')
        if len(p) < 50 or random.random() > 0.7:
            continue
        ko = random.choice(kos[:4000]) if random.random() < 0.8 else random.choice(kos)
        for copy in range(random.choice([1, 1, 2])):
            rate = random.choice([0.03, 0.1, 0.25, 0.4, 0.6])
            m = ''.join(random.choice(AA) if random.random() < rate else ch for ch in p)
            # occasional length change (coverage filter, length-ratio seeding)
            if random.random() < 0.1:
                m = m[: int(len(m) * random.choice([0.4, 0.7]))]
            n += 1
            refs.append(('UniRef100_%07d~%s' % (n, ko), m))
for i in range(3000):  # decoys
    n += 1
    refs.append(('UniRef100_%07d~%s' % (n, random.choice(kos)), ''.join(random.choice(AA) for _ in range(random.randint(80, 600)))))
os.makedirs(out + '/db/checkm2/CheckM2_database', exist_ok=True)
write(out + '/ref.faa', refs)
subprocess.run(['diamond', 'makedb', '--in', out + '/ref.faa', '-d', out + '/db/checkm2/CheckM2_database/uniref100.KO.1', '--quiet'], check=True)
print('references', len(refs))
