"""Compare CheckM2's feature vectors (--dbg_vectors pickle) and raw model
outputs with Sentinel's dumps (SENTINEL_CHECKM2_DEBUG)."""
import sys, os, warnings, pickle
os.environ.setdefault('KERAS_BACKEND', 'jax')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
sys.path.insert(0, '/home/claude/refsrc/checkm2')
from checkm2 import modelProcessing, modelPostprocessing
pkl, sdir = sys.argv[1], sys.argv[2]
fv = pd.read_pickle(pkl)
names = fv['Name'].values
X = fv.iloc[:, 1:].values.astype(float)
mp = modelProcessing.modelProcessor(1)
gen, cont = mp.run_prediction_general(X)
spec, scaled = mp.run_prediction_specific(X, 20021)
post = modelPostprocessing.modelProcessor(1)
_, _, chosen, cos = post.calculate_general_specific_ratio(X[:, 20], scaled, gen, cont, spec)
for i, n in enumerate(names):
    sv = np.array([float(l) for l in open(os.path.join(sdir, n + '.features.txt'))])
    nd = int(np.sum(sv != X[i]))
    pred = dict(l.split('\t') for l in open(os.path.join(sdir, n + '.prediction.txt')).read().split('\n') if l)
    print(f"{n:8s} features differing: {nd}/{len(sv)}  "
          f"general {gen[i]!r} vs {pred['general']}  cont {cont[i]!r} vs {pred['contamination']}  "
          f"specific {float(spec[i])!r} vs {pred['specific']}  cosine {cos[i]!r} vs {pred['cosine']}")
