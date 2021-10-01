import sys
import vtkcommon
import vtk
import h5py
import pandas as pd
import numpy as np


# create an average chemical profile for each region


def main():
    salm_density = pd.read_csv('data/bacterial-density.csv')
    # line of best fit for the overall bacterial density
    m, b = np.polyfit(salm_density['time'], salm_density['density'], 1)

    region_accs = [
        pd.read_csv('data/relative-body.csv'),
        pd.read_csv('data/relative-core.csv'),
        pd.read_csv('data/relative-transition.csv'),
        pd.read_csv('data/relative-edge.csv'),
    ]

    region_densities = [[] for _ in range(len(region_accs))]
    for i in range(len(region_accs[0])):
        time = region_accs[0].iloc[i]['time']
        relative_total = sum([x.iloc[i]['relative-acc']  for x in region_accs])
        bacterial_density = time * m + b
        for j,r in enumerate(region_accs):
            region_density_frac = (r.iloc[i]['relative-acc'] / relative_total) \
                              #* bacterial_density
            region_densities[j].append(region_density_frac)
    for i,r in enumerate(region_accs):
        r['density_frac'] = region_densities[i]
        print(r)

    
    timesteps = [0,12,48,96]
    # generate vtk files for timesteps 0, 12, 48, 96
    with h5py.File(sys.argv[1], 'r') as tumor_h5:
        # Edge
        regions = np.array(tumor_h5['out0540']['cells']['region'])

        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(sys.argv[2])
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        polydata = reader.GetOutput()
        num_cells = polydata.GetNumberOfPoints()

        for i,t in enumerate(timesteps):
            # for t=0 just put density = 0 for all cells
            data = np.zeros((num_cells,1))
            if t != 0:
                for j, r in enumerate(region_accs):
                    data[regions == j] = r.iloc[i-1]['density_frac']
            polydata.GetPointData().AddArray(vtkcommon.asVtkArray(data, 'salm', vtk.vtkFloatArray))

            writer = vtk.vtkPolyDataWriter()
            writer.SetInputData(polydata)
            writer.SetFileName(f'tumor-{t}.vtk')
            writer.Write()



if __name__=='__main__':
    sys.exit(main())
