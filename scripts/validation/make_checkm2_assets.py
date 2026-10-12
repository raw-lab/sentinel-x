"""Convert CheckM2 1.1.0's specific model, scaler and reference data into
Sentinel-X's compact embedded assets.

  specific_model.bin.z : 4x(Conv1D W,b; BN gamma,beta,mean,var), 2x Dense W,b,
                         then the MinMaxScaler scale_ and min_ (all 21,241)
  ref_rsdata.bin.z     : min_ref_rsdata_v1.npz (first 20,021 columns) as
                         integer feature counts: metadata block, column-major
                         KO presence bitmap, KO counts
Values are reconstructed exactly as count*scale_ + min_ (verified here).
"""
import sys, zlib, struct, pickle, zipfile, io, warnings
import numpy as np, scipy.sparse as sp, h5py
warnings.filterwarnings("ignore")
src, out = sys.argv[1], sys.argv[2]
M = src + '/models/'; D = src + '/data/'
sc = pickle.load(open(M + 'scaler.sav', 'rb'))

# ---- specific model + scaler
z = zipfile.ZipFile(M + 'specific_model_COMP.keras')
h = h5py.File(io.BytesIO(z.read('model.weights.h5')), 'r')
L = h['layers']
buf = bytearray(b'SNN1')
def put(a, dt):
    a = np.ascontiguousarray(a, dtype=dt)
    buf.extend(struct.pack('<I', a.size)); buf.extend(a.tobytes())
for conv, bn in [('conv1d', 'batch_normalization'), ('conv1d_1', 'batch_normalization_1'),
                 ('conv1d_2', 'batch_normalization_2'), ('conv1d_3', 'batch_normalization_3')]:
    W = L[conv]['vars']['0'][()]; b = L[conv]['vars']['1'][()]
    buf.extend(struct.pack('<III', *W.shape)); put(W, '<f4'); put(b, '<f4')
    for k in range(4):
        put(L[bn]['vars'][str(k)][()], '<f4')
for d in ['dense', 'dense_1']:
    W = L[d]['vars']['0'][()]; b = L[d]['vars']['1'][()]
    buf.extend(struct.pack('<II', *W.shape)); put(W, '<f4'); put(b, '<f4')
put(sc.scale_, '<f8'); put(sc.min_, '<f8')
open(out + '/specific_model.bin.z', 'wb').write(zlib.compress(bytes(buf), 9))

# ---- reference data
NC, NM = 20021, 22
r = sp.load_npz(D + 'min_ref_rsdata_v1.npz')[:, :NC].tocsr(); r.sort_indices()
scale = sc.scale_[:NC]; mn = sc.min_[:NC]
cnt = np.rint((r.data - mn[r.indices]) / scale[r.indices]).astype(np.int64)
assert np.all(cnt * scale[r.indices] + mn[r.indices] == r.data), 'values are not count*scale+min'
assert cnt.min() >= 1
rows = r.shape[0]
def varint(n, o):
    while True:
        b = n & 0x7f; n >>= 7
        if n: o.append(b | 0x80)
        else: o.append(b); break
o = bytearray(b'SREF1')
o.extend(struct.pack('<III', rows, NC, NM))
dense_meta = np.zeros((rows, NM), dtype=np.int64)
bm = np.zeros((rows, NC - NM), dtype=np.uint8)
kocounts = []
for i in range(rows):
    s, e = r.indptr[i], r.indptr[i + 1]
    for jj in range(s, e):
        j = r.indices[jj]
        if j < NM:
            dense_meta[i, j] = cnt[jj]
        else:
            bm[i, j - NM] = 1
            kocounts.append(cnt[jj])
for i in range(rows):
    for j in range(NM):
        varint(int(dense_meta[i, j]), o)
o.extend(np.packbits(bm.T, axis=1).tobytes())
for c in kocounts:
    varint(int(c) - 1, o)
open(out + '/ref_rsdata.bin.z', 'wb').write(zlib.compress(bytes(o), 9))

# metadata zero entries are absent from the CSR: check explicitly
nz_meta = sum(1 for i in range(rows) for j in range(NM) if dense_meta[i, j] == 0)
print('rows', rows, 'nnz', r.nnz, 'zero metadata cells', nz_meta)
