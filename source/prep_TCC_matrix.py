import os
import sys, gc

json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print("ERROR: Please provide path to a valid config.json file...")
    print(sys.argv[1])
    exit(1)
    
import json
with open(json_path) as json_file:
    parameter = json.load(json_file)

print("Number of threads", parameter["NUM_THREADS"])

# Load dataset



import numpy as np
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize

#matrix.ec file
ecfile_dir = parameter["kallisto"]["TCC_output"]+"matrix.ec"
tsvfile_dir = parameter["kallisto"]["TCC_output"]+"matrix.tsv"

print("Loading TCCs..")

COOinput = np.loadtxt( tsvfile_dir, delimiter='\t' , dtype=int)
rows,cols,data = COOinput.T
nonzero_ec = np.unique(rows)
map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
TCCmatrix   = coo_matrix( (data.astype(float),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) ) 

NUM_OF_CELLS = TCCmatrix.shape[1]
print("NUM_OF_CELLS =", NUM_OF_CELLS)
      
T = TCCmatrix.tocsr()
T_norm = normalize(T, norm='l1', axis=0) 
T_normT = T_norm.transpose()
del TCCmatrix;
_ = gc.collect()



# Pairwise_distances


from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import *
from scipy.stats import entropy

def L1_distance(p,q):
    return cityblock(p,q).sum()

# def jensen_shannon(p, q):
#     m=0.5*p+0.5*q
#     p = np.transpose(p[p > 0])
#     q = np.transpose(q[q > 0])
#     m = np.transpose(m[m > 0])
#     return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))

num_of_threads = parameter["NUM_THREADS"]
print("Calculating pairwise L1 distances... ( num_threads =",num_of_threads,")")

# D_js = pairwise_distances(T_normT,metric=jensen_shannon,n_jobs=num_of_threads)
D_l1 = pairwise_distances(T_normT,metric=L1_distance,n_jobs=num_of_threads)

print("writing data...")

#Save data
import pickle

with open(parameter["SAVE_DIR"]+"TCC_matrix.dat", 'wb') as f:
    pickle.dump(T,f)
with open(parameter["SAVE_DIR"]+"pwise_dist_L1.dat", 'wb') as f:
    pickle.dump(D_l1,f)
with open(parameter["SAVE_DIR"]+"nonzero_ec.dat", 'wb') as f:
    pickle.dump(nonzero_ec,f)

print("DONE.")





