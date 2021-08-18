import h5py
import sys
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import RobustScaler



variables = [
    'cells/cell_phase',
    'cells/cell_phase_age',
    'cells/concAcL_ex',
    'cells/pH_ex',
    'cells/glucose_ex',
    'cells/cell_o2_consumption_rate',
    'cells/cell_radii',
    'cells/distance_to_nearest_vessel',
    'vbl/isonAS',
    'vbl/age',
    'vbl/age_mother',
    'vbl/volume_extra'
]


def main():
    with h5py.File(sys.argv[1], 'r') as tumor_h5:
        num_cells = int(tumor_h5['out0540/vbl'].attrs['ncells'])
        num_variables = len(variables)

        # load the data into np array
        X = np.empty((num_cells, num_variables))
        for i,v in enumerate(variables):
            assert tumor_h5['out0540'][v].shape == (num_cells, 1)
            X[:, i] = np.reshape(tumor_h5['out0540'][v], (num_cells,))

        # scale the data
        X = RobustScaler().fit_transform(X)

        # calculate distance from (0,0,0) for each cell
        y = np.linalg.norm(tumor_h5['out0540/cells/cell_center_pos'], axis=1)

        # fit linear regression and determine variables most predictive of location
        linreg = LinearRegression().fit(X, y)
        var_importance = sorted(zip(range(num_variables), np.abs(linreg.coef_)),
                                key=lambda x:x[1],
                                reverse=True)

        for x in var_importance:
            print(f"{variables[x[0]]:32}: {x[1]}")




if __name__=='__main__':
    sys.exit(main())


