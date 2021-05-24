import os
import shutil
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import imageio


def main():
    file = sys.argv[1]          # path to h5 file
    timestep = sys.argv[2]      # e.g. out0002

    h5file = h5py.File(file)

    datasets =  [
        h5file[timestep]['oxy_source_lin'],
        h5file[timestep]['vessel_volume_fraction'],
        h5file[timestep]['fieldGf'],
        h5file[timestep]['fieldOxy'],
        h5file[timestep]['tumor']['conc'],
        h5file[timestep]['tumor']['necro'],
        h5file[timestep]['tumor']['ls'],
        h5file[timestep]['tumor']['ptc'],
        h5file[timestep]['tumor']['sources'],
        h5file[timestep]['tumor']['press'],
    ]

    shape = (74,74,64)

    for ds in datasets:
        assert(ds.shape == shape)

    vranges = [(np.min(ds), np.max(ds)) for ds in datasets]

    shutil.rmtree('frames', ignore_errors=True)
    os.makedirs('frames')

    for i in range(shape[2]):
        fig = plt.figure(figsize=(20, 10), dpi=100)
        plt.title(f'z={i}\n\n')
        plt.axis("off")
        for j,ds in enumerate(datasets):
            f = fig.add_subplot(2, 5, j+1)
            f.axis('off')
            f.title.set_text(ds.name.split('/')[-1])
            im1 = f.imshow(ds[:,:,i], vmin=vranges[j][0], vmax=vranges[j][1])
            cb = fig.colorbar(im1, location='bottom')
            cb.ax.tick_params(labelsize=8)
        plt.savefig(f'frames/{i}.png')
        plt.close()

    with imageio.get_writer(f'{timestep}.gif', mode='I', duration=0.25) as writer:
        for i in range(shape[2]):
            image = imageio.imread(f'frames/{i}.png')
            writer.append_data(image)

if __name__ == '__main__':
    sys.exit(main())

