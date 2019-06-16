import struct
import numpy as np
import h5py
import os


def getParticles(fname, hdf5 = True):
    if (hdf5):
        # hd5 file
        with h5py.File(fname, 'r') as file:
            keys = list(file.keys())
            species = np.unique([int(s) for s in ''.join(keys) if s.isdigit()])
            nspec = len(species)
            variables = np.unique([''.join([i for i in k if not i.isdigit()]) for k in keys])
            nvars = len(variables)
            data = {}
            for s in range(nspec):
                data[str(s + 1)] = {}
                for i in range(nvars):
                    (data[str(s + 1)])[variables[i]] = file[variables[i] + str(s + 1)][:]
    else:
        # !--- PRTL.TOT.***** structure ----------------------------------!
        # ! HEADER:                                                      _
        # !   timestep......................[4 bytes]                     |
        # !   # of cpus.....................[4 bytes]                     |
        # !   # of species..................[4 bytes]                     |
        # !   # of variables................[4 bytes]                     |- disp_header
        # !   variable names................[#var * 5 bytes]              |
        # !   variable types................[#var * 5 bytes]              |
        # !   # of particles per species....                              |
        # !   ....summed over all ranks.....[#spec * 4 bytes]            _|
        # ! BODY:
        # !   species = 1...............[#of parts of species=1 * # of variables * 4 bytes]
        # !     var = 1.................[#of parts per species * 4 bytes]
        # !       rank = 1..............[#of parts on rank=1 per species * 4 bytes]
        # !         XXX.................[4 bytes]
        # !         XXX.................[4 bytes]
        # !         ....................
        # !         XXX.................[4 bytes]
        # !       rank = 2..............[#of parts on rank=1 per species * 4 bytes]
        # !       ......................
        # !       rank = N..............[#of parts on rank=N per species * 4 bytes]
        # !     var = 2.................[#of parts per species * 4 bytes]
        # !     ........................
        # !     var = V.................[#of parts per species * 4 bytes]
        # !   species = 2...............[#of parts of species=2 * # of variables * 4 bytes]
        # !   ..........................
        # !   species = S...............[#of parts of species=S * # of variables * 4 bytes]
        # !   ..........................
        # !...............................................................!
        # binary file
        with open(fname, mode='rb') as file:
            fileContent = file.read()
            data = {}
            header_start = 4
            timestep, mpi_size, nspec, nvars = struct.unpack("i" * header_start, fileContent[:4 * header_start])
            data['timestep'] = timestep
            variables = []
            variable_types = []
            read_ptr = 4 * header_start
            for i in range(nvars):
                varn = struct.unpack("s"*5, fileContent[read_ptr : read_ptr + 5])
                varn = b''.join(varn).replace(b' ', b'').decode('ascii')
                read_ptr += 5
                variables.append(varn)
            for i in range(nvars):
                vart = struct.unpack("s"*5, fileContent[read_ptr : read_ptr + 5])
                vart = b''.join(vart).replace(b' ', b'').decode('ascii')
                read_ptr += 5
                if (vart == 'real'):
                    variable_types.append("f")
                elif (vart == 'int'):
                    variable_types.append("i")
            nprt = struct.unpack("i" * nspec, fileContent[read_ptr : read_ptr + nspec * 4])
            read_ptr += nspec * 4
            for s in range(nspec):
                data[str(s + 1)] = {}
                for i in range(nvars):
                    (data[str(s + 1)])[variables[i]] = np.array(struct.unpack(variable_types[i] * nprt[s], fileContent[read_ptr : read_ptr + nprt[s] * 4]))
                    read_ptr += nprt[s] * 4
    return data

def convertPartsToHdf5(fname, delete_original = True):
    particles = getParticles(fname)
    if delete_original:
        newname = fname
        os.remove(fname)
    else:
        newname = fname.replace("tot", "hdf5")
    nspec = len(particles) - 1
    varr = list(particles['0'].keys())
    hf = h5py.File(newname, 'w')
    hf.create_dataset('timestep', data=particles['timestep'])
    for i in range(nspec):
        group = hf.create_group(str(i))
        for v in varr:
            dset = group.create_dataset(v, data=particles[str(i)][v])
    hf.close()

def getFields(fname, hdf5 = True, nodes = False):
    if hdf5:
        # hdf5 file
        with h5py.File(fname, 'r') as file:
            keys = list(file.keys())
            data = {}
            for key in keys:
                data[key] = file[key][:]
    else:
        # !--- FLDS.TOT.***** structure ----------------------------------!
        # ! HEADER:                                                      _
        # !   timestep......................[4 bytes]                     |
        # !   # of cpus.....................[4 bytes]                     |
        # !   # of fields...................[4 bytes]                     |
        # !   dimensions..[sx,sy,sz]........[3 * 4 bytes]                 |- disp_header
        # !   field names...................[#flds * 5 bytes]             |
        # !   meshblock dimensions..........[# of cpus * 6 * 4 bytes]    _|
        # ! BODY:
        # !   field = 1.................[fx * fy * fz * 4 bytes]
        # !     rank = 1................[fx * fy * fz (for rank = 1) * 4 bytes]
        # !       XXX...................[4 bytes]
        # !       XXX...................[4 bytes]
        # !       ......................
        # !       XXX...................[4 bytes]
        # !     rank = 2................[fx * fy * fz (for rank = 2) * 4 bytes]
        # !     ........................
        # !     rank = N................[fx * fy * fz (for rank = N) * 4 bytes]
        # !   field = 1.................[fx * fy * fz * 4 bytes]
        # !   ..........................
        # !   field = F.................[fx * fy * fz * 4 bytes]
        # !   ..........................
        # !...............................................................!
        # sx,sy,sz -> original dimensions
        # fx,fy,fz -> downsampled dimensions
        # binary file
        def getGlobalS(x0, siz):
            siz = np.array([x for _,x in sorted(zip(x0,siz))])
            x0 = np.array(sorted(x0))
            mmin = siz[0]
            curr = x0[0]
            for rnk in range(4):
                if x0[rnk] > curr:
                    mmin += siz[rnk]
                    curr = x0[rnk]
            return mmin
        with open(fname, mode='rb') as file:
            fileContent = file.read()
            data = {}
            header_start = 3 + 3
            timestep, mpi_size, nflds, dimx, dimy, dimz = struct.unpack("i" * header_start, fileContent[:4 * header_start])
            data['timestep'] = timestep
            variables = []
            read_ptr = 4 * header_start
            for i in range(nflds):
                varn = struct.unpack("s"*5, fileContent[read_ptr : read_ptr + 5])
                varn = b''.join(varn).replace(b' ', b'').decode('ascii')
                read_ptr += 5
                variables.append(varn)
            x_y_z_list = []
            sx_sy_sz_list = []
            for i in range(mpi_size):
                x_y_z = struct.unpack("i"*3, fileContent[read_ptr : read_ptr + 4 * 3])
                read_ptr += 4 * 3
                sx_sy_sz = struct.unpack("i"*3, fileContent[read_ptr : read_ptr + 4 * 3])
                read_ptr += 4 * 3
                x_y_z_list.append(x_y_z); sx_sy_sz_list.append(sx_sy_sz)
            data['mblocks_xyz'] = np.array(x_y_z_list)
            data['mblocks_sxyz'] = np.array(sx_sy_sz_list)
            # getting global domain sizes from the local meshblocks
            sx_glob = getGlobalS(np.array(x_y_z_list)[:,0], np.array(sx_sy_sz_list)[:,0])
            sy_glob = getGlobalS(np.array(x_y_z_list)[:,1], np.array(sx_sy_sz_list)[:,1])
            sz_glob = getGlobalS(np.array(x_y_z_list)[:,2], np.array(sx_sy_sz_list)[:,2])
            # saving fields
            for f in range(nflds):
                fld_glob = np.zeros((sx_glob, sy_glob, sz_glob))
                for rnk in range(mpi_size):
                    x0, y0, z0 = x_y_z_list[rnk]
                    sx, sy, sz = sx_sy_sz_list[rnk]
                    fld = np.array(struct.unpack("f" * sx * sy * sz, fileContent[read_ptr : read_ptr + 4 * sx * sy * sz]))
                    read_ptr += 4 * sx * sy * sz
                    fld = fld.reshape(sx, sy, sz) # then `fld[xi,yi,zi]` is the field at `xi,yi,zi`
                    fld_glob[x0:x0+sx,y0:y0+sy,z0:z0+sz] = fld
                data[variables[f]] = np.array(fld_glob)
            # saving grid as a field
            fld_glob = np.zeros((sx_glob, sy_glob, sz_glob, 3))
            for rnk in range(mpi_size):
                x0, y0, z0 = x_y_z_list[rnk]
                sx, sy, sz = sx_sy_sz_list[rnk]
                xyz_grid = [[[(i + x0, j + y0, k + z0) for k in range(sz)] for j in range(sy)] for i in range(sx)]
                fld_glob[x0:x0+sx,y0:y0+sy,z0:z0+sz] = xyz_grid
            if (nodes):
                x_ = fld_glob[:,0,0,0]
                x_ = np.append(x_, x_[-1] + (x_[-1] - x_[-2]))
                y_ = fld_glob[0,:,0,1]
                y_ = np.append(y_, y_[-1] + (y_[-1] - y_[-2]))
                z_ = fld_glob[0,0,:,2]
                if len(z_) > 1:
                    z_ = np.append(z_, z_[-1] + (z_[-1] - z_[-2]))
                else:
                    z_ = np.append(z_, z_[-1] + 1)
                data['x'] = x_
                data['y'] = y_
                data['z'] = z_
            else:
                data['xyz'] = np.array(fld_glob)
            if (len(fileContent) != read_ptr):
                print ("WRONG reading!")
    return data
# usage example for 2D uniform grid:
# ```
#   field_data = isolde.getFields("flds.tot.00000", True)
#   x_ = field_data['x']
#   y_ = field_data['y']
#   x_, y_ = np.mgrid[x_[0]: x_[-1] + 1 : x_[1] - x_[0],
#                   y_[0]: y_[-1] + 1 : y_[1] - y_[0]]
#   ex_ = field_data['ex'][:,:,0]
#   plt.pcolor(x_, y_, ex_) # <- 2D plot
# ```

def getSpectra(fname, hdf5 = True):
    if hdf5:
        # hdf5 output
        with h5py.File(fname, 'r') as file:
            keys = list(file.keys())
            species = np.unique([int(s) for s in ''.join(keys) if s.isdigit()])
            nspec = len(species)
            data = {}
            for s in range(nspec):
                data[str(s + 1)] = {}
                (data[str(s + 1)])['bn'] = np.exp(file['e' + str(s + 1)][:])
                (data[str(s + 1)])['cnt'] = file['n' + str(s + 1)][:]
    else:
        with open(fname, mode='rb') as file:
            # !--- SPEC.TOT.***** structure ----------------------------------!
            # ! HEADER:
            # !   timestep......................[4 bytes]
            # !   # of species..................[4 bytes]
            # !   # of bins.....................[4 bytes]
            # !   MIN energy....................[4 bytes]
            # !   MAX energy....................[4 bytes]
            # ! BODY:
            # !   species = 1...............[# of bins * 4 bytes]
            # !     bin = 1.................[4 bytes]
            # !     bin = 2.................[4 bytes]
            # !     ........................
            # !     rank = N................[4 bytes]
            # !   species = 2...............[# of bins * 4 bytes]
            # !   ..........................
            # !   species = S...............[# of bins * 4 bytes]
            # !   ..........................
            # !...............................................................!
            # binary output
            fileContent = file.read()
            data = {}
            timestep, nspec, nbins = struct.unpack("i" * 3, fileContent[:4 * 3])
            read_ptr = 4 * 3
            emin, emax = struct.unpack("f" * 2, fileContent[read_ptr : read_ptr + 4 * 2])
            read_ptr += 4 * 2
            data['timestep'] = timestep
            data['nspec'] = nspec
            emin = np.round(np.log10(np.exp(emin)), 4)
            emax = np.round(np.log10(np.exp(emax)), 4)
            data['bins'] = np.logspace(emin, emax, nbins)
            for s in range(nspec):
                data['spec' + str(s + 1)] = np.array(struct.unpack("i" * nbins, fileContent[read_ptr : read_ptr + 4 * nbins]))
                read_ptr += 4 * nbins
    return data

def getDomains(fname):
    with h5py.File(fname, 'r') as file:
        data = {}
        for k in file.keys():
            data[k] = file[k][:]
    return data

# easy plotting functions
def plot2DField(ax, x, y, field,
                title='field', cmap='jet',
                vmin=None, vmax=None,
                typ='lin', **kwargs):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    sx, sy = field.shape
    if not vmin:
        vmin = field.min()
    if not vmax:
        vmax = field.max()
    if typ == 'lin':
        im = ax.pcolormesh(x, y, field, cmap=cmap, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    elif typ == 'log':
        im = ax.pcolormesh(x, y, field, cmap=cmap, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
    elif typ == 'sym':
        vmax = max(np.abs(vmin), vmax)
        im = ax.pcolormesh(x, y, field, cmap=cmap,
                           norm=mpl.colors.SymLogNorm(vmin=-vmax, vmax=vmax,
                                                      linthresh=kwargs['lth'],
                                                      linscale=kwargs['lsc']))
    ax.set_aspect(1)
    ax.set_xlim(0, sx - 1)
    ax.set_ylim(0, sy - 1)
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    else:
        ax.set_xlabel('x')
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])
    else:
        ax.set_ylabel('y')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    ax.set_title(title)
    plt.colorbar(im, cax=cax)

def plot2DScatterParticles(ax, sx, sy, x_list, y_list,
                           label='particles', legend=True,
                           color='black', **kwargs):
    ax.scatter(x_list, y_list, c=color, label=label)
    ax.set_aspect(1)
    if legend:
        ax.legend()
    ax.set_xlim(0, sx - 1)
    ax.set_ylim(0, sy - 1)

def plot2DDomains(ax, domain_data,
                  color='red', **kwargs):
    from matplotlib.patches import Rectangle
    x0_list = domain_data['x0']
    y0_list = domain_data['y0']
    sx_list = domain_data['sx']
    sy_list = domain_data['sy']
    for x0, y0, sx, sy in zip(x0_list, y0_list, sx_list, sy_list):
        rect = Rectangle((x0, y0), sx, sy,
                         edgecolor=color, facecolor='none')
        ax.add_patch(rect)
