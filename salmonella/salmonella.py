import h5py
import sys
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import RobustScaler
import skfuzzy
import math

# using this paper as reference: https://www.nature.com/articles/cgt201110

# Region definitions:
# - Edge
# - Body
# - Transition
# - Core

# Salmonella concentration will depend on:
# - Region
# - Time

# Aim: fuzzy classify each cell to regions

# Steps:
# 1) look at all of the potentially relevant variables
# 3) calculate distance to center of tumor for each cell
# 4) fit a linear regression model (normalized) that relates those variables to the distance from the centre of the tumor
# 5) Pick the top N variables by their regression coefficients
# 6) perform fuzzy-c-means clustering with c = 4 to assign a region score to each cell

# We'll define 2 distinct regions 

variables = [
    'cells/concAcL_ex',
    'cells/glucose_ex',
]


def main():
    with h5py.File(sys.argv[1], 'r') as tumor_h5:
        num_cells = int(tumor_h5['out0540/vbl'].attrs['ncells'])
        num_variables = len(variables)
        X = np.empty((num_cells, num_variables+1))

        # load the data into an np array
        for i,v in enumerate(variables):
            assert tumor_h5['out0540'][v].shape == (num_cells, 1)
            X[:, i] = np.reshape(tumor_h5['out0540'][v], (num_cells,))
        
        # include distance from centre as one of the clustering parameters
        X[:, num_variables] = np.linalg.norm(tumor_h5['out0540/cells/cell_center_pos'], axis=1)

        X = RobustScaler().fit_transform(X)

        cntr, u, u0, d, jm, p, fpc = skfuzzy.cluster.cmeans(X.T, 2, 2, 
                                                             error=0.005,
                                                             maxiter=1000,
                                                             init=None)
        scores = u.T
        regions = np.argmax(scores, axis=1)

        # Transition
        threshold = 0.1
        for i,(a,b) in enumerate(scores):
            if abs(a-b) < threshold:
                regions[i] = 2

        l_regions = np.loadtxt('regions.txt', dtype=int)
        regions[np.argwhere(l_regions == 3)] = 3

        # Edge
#        # Create array of all cells (their positions) that are on the outside surface (isonAS)
#        edge_indices = np.nonzero(np.reshape(tumor_h5['out0540/vbl/isonAS'], (num_cells,)))[0]
#        edge_pos = tumor_h5['out0540']['cells']['cell_center_pos'][edge_indices]
#
#        cell_pos = tumor_h5['out0540/cells/cell_center_pos']
#
#        def closest_node_idx(node, nodes):
#            nodes = np.asarray(nodes)
#            deltas = nodes - node
#            dist_2 = np.einsum('ij,ij->i', deltas, deltas)
#            return np.argmin(dist_2)
#
#        count = 0
#        threshold = 20
#        for i,cp in enumerate(cell_pos):
#            closest = edge_pos[closest_node_idx(cp, edge_pos)]
#            distance = np.linalg.norm(cp-closest)
#            if distance < threshold:
#                count += 1
#                regions[i] = 3
#
#        print(count / num_cells)
#        np.savetxt('regions.txt', regions, fmt='%d')

        with h5py.File("tumor-regions.h5", 'w') as new_h5:
            for a in tumor_h5.attrs:
                new_h5.attrs[a] = tumor_h5.attrs[a]
            new_h5.copy(tumor_h5['out0540'], new_h5)
            new_h5.copy(tumor_h5['field_ld'], new_h5)
            new_h5.copy(tumor_h5['last_state'], new_h5)
            new_h5.copy(tumor_h5['parameters'], new_h5)

            region_ds = new_h5['out0540']['cells'].create_dataset('region',
                                                               (num_cells, 1),
                                                               dtype='i8',
                                                               data=regions)
            
            
  

if __name__=='__main__':
    sys.exit(main())

