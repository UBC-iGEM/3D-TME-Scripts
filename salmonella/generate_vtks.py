import sys
import vtk
import vtkcommon
import numpy as np


def main():
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(sys.argv[1])
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    polydata = reader.GetOutput()

    test = np.zeros((421646,))
    polydata.GetPointData().AddArray(vtkcommon.asVtkArray(test, 'test', vtk.vtkFloatArray))

    print(polydata)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName('test')
    writer.Write()




if __name__=='__main__':
    sys.exit(main())
