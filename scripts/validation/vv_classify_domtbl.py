"""Classify contigs from a --domtblout with viralVerify's own functions."""
import sys, os, types, collections
vv = '/home/claude/refsrc/viralVerify-master/viralVerify-master'
sys.path.insert(0, vv + '/scripts')
src = open(vv + '/bin/viralverify').read().replace('__file__', repr(vv + '/bin/viralverify'))
mod = types.ModuleType('vv'); mod.__dict__['__name__'] = 'vv'
exec(compile(src.split("if __name__")[0], 'viralverify', 'exec'), mod.__dict__)
from check_circular import check_circular
domtbl, fasta = sys.argv[1], sys.argv[2]
circ = check_circular(fasta, '/tmp/vvcirc')
ft = [i.strip().split(' ', 1) for i in mod.get_table_from_tblout(domtbl)]
k = mod.naive_bayes([i[1] for i in ft], 7)
res = {ft[i][0]: [k[i][0], k[i][3], ft[i][1]] for i in range(len(k))}
for cid, (L, c) in circ.items():
    if cid in res:
        p = res[cid][0]
        if p == 'Uncertain - too short' and (L > 3000 or c == '+'): p = 'Uncertain - viral or bacterial'
        print(cid, p, res[cid][1], ' '.join(sorted(res[cid][2].split())), sep='\t')
    else:
        print(cid, 'Uncertain - viral or bacterial' if (L > 3000 or c == '+') else 'Uncertain - too short', '-', '', sep='\t')
