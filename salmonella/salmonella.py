import h5py
import sys
import numpy as np
import skfuzzy
from sklearn.preprocessing import RobustScaler

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
    'cells/pH_ex',
]


EDGE_SIZE = 0.05

# let's do k means with 2 regions.
# we'll define the transition zone as all cells that are close to both

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

        cntr, u, u0, d, jm, p, fpc = skfuzzy.cluster.cmeans(X.T, 4, 2, 
                                                            error=0.005,
                                                            maxiter=1000,
                                                            init=None)
        scores = u.T
        regions = np.argmax(scores, axis=1)

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

