"""Local mirror of the database sources for downloader tests."""
import os, sys, json, hashlib, tarfile, gzip, io, shutil
site = sys.argv[1]
W = '/home/claude/work'
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()
def sha256(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
# CheckV
open(site + '/checkv/CURRENT_RELEASE.txt', 'w').write('checkv-db-v1.5-mini\n')
with tarfile.open(site + '/checkv/checkv-db-v1.5-mini.tar.gz', 'w:gz') as t:
    t.add(W + '/h2h_checkv/minidb/checkv-db-v1.5-mini', arcname='checkv-db-v1.5-mini')
# CheckM2 backpack
dmnd = W + '/h2h_checkm2/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd'
contents = {'version': '3', 'zenodo_backpack_version': 1, 'CheckM2_database/uniref100.KO.1.dmnd': md5(dmnd)}
cj = json.dumps(contents).encode()
tgz = site + '/files/checkm2_database.tar.gz'
with tarfile.open(tgz, 'w:gz') as t:
    ti = tarfile.TarInfo('CONTENTS.json'); ti.size = len(cj); t.addfile(ti, io.BytesIO(cj))
    t.add(dmnd, arcname='CheckM2_database/uniref100.KO.1.dmnd')
port = sys.argv[2]
rec = {'metadata': {'version': '3'}, 'files': [{'key': 'checkm2_database.tar.gz',
       'links': {'self': f'http://127.0.0.1:{port}/files/checkm2_database.tar.gz'}, 'checksum': 'md5:' + md5(tgz)}]}
open(site + '/api/records/14897628', 'w').write(json.dumps(rec))
# Pfam
hmm = W + '/vv_h2h/mini_pfam.hmm'
with gzip.open(site + '/pfam/Pfam-A.hmm.gz', 'wb') as g: g.write(open(hmm, 'rb').read())
with gzip.open(site + '/pfam/Pfam.version.gz', 'wb') as g: g.write(b'Pfam release       : 99.0\nPfam-A families    : 548\n')
open(site + '/pfam/md5_checksums', 'w').write(md5(site + '/pfam/Pfam-A.hmm.gz') + '  Pfam-A.hmm.gz\n' + '0' * 32 + '  Pfam-B.fasta.gz\n')
# BUSCO
tb = site + '/busco/lineages/mini_odb10.2026-10-08.tar.gz'
with tarfile.open(tb, 'w:gz') as t:
    t.add(W + '/busco_h2h/lineages/mini_odb10', arcname='mini_odb10')
open(site + '/busco/file_versions.tsv', 'w').write(
    f'mini_odb10\t2026-10-08\t{md5(tb)}\tProkaryota\tlineages\nbacteria_odb12.2\t2026-05-22\te21960d275e09687a86aace486f2886a\tProkaryota\tlineages\n')
print('dmnd sha256', sha256(dmnd))
