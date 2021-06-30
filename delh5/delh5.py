import h5py
import sys




def main():
    with h5py.File(sys.argv[1], 'r') as orig:
        with h5py.File("out.h5", 'w') as new:
            for a in orig.attrs:
                new.attrs[a] = orig.attrs[a]
            new.copy(orig['out0540'], new)
            new.copy(orig['field_ld'], new)
            new.copy(orig['last_state'], new)
            new.copy(orig['parameters'], new)







if __name__=='__main__':
    main()

