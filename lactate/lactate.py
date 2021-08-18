import h5py
import sys

# https://github.com/thierry3000/VBL/blob/master/src/CellsSystem.cpp#L708

BufCapEnv = 0.584e-3
# 89.07 g/mol is the MW of lactate
# returns extracellular [lactate] in mM
 
def get_AcL(pH_ex):
    return (7.5443 - pH_ex) * BufCapEnv * 1.0e6 * (1/89.07)

def main():
    with h5py.File(sys.argv[1], 'r') as orig:
        with h5py.File("tumor-lac.h5", 'w') as new:
            for a in orig.attrs:
                new.attrs[a] = orig.attrs[a]
            new.copy(orig['out0540'], new)
            new.copy(orig['field_ld'], new)
            new.copy(orig['last_state'], new)
            new.copy(orig['parameters'], new)
            cells_group = new['out0540']['cells']
            pH_ex = cells_group['pH_ex']
            cells_group.create_dataset('AcL_ex', data=pH_ex)
            AcL_ex = cells_group['AcL_ex']
            AcL_ex[:] = [get_AcL(pH_ex[i]) for i in range(AcL_ex.shape[0])]

if __name__=='__main__':
    sys.exit(main())

