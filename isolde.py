import struct
import numpy as np
import h5py
import os

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
def getParticles(fname):
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
def getFields(fname):
    def globalizeFld(fld_lst, tuples=False):
        if tuples:
            fld0 = np.zeros((dimx, dimy, dimz, 3))
        else:
            fld0 = np.zeros((dimx, dimy, dimz))
        def findRank(i, j, k):
            for rnk in range(mpi_size):
                x0, y0, z0 = x_y_z_list[rnk]
                sx, sy, sz = sx_sy_sz_list[rnk]
                if ((i >= x0) and (i < x0 + sx) and (j >= y0) and (j < y0 + sy) and (k >= z0) and (k < z0 + sz)):
                    return rnk
            return -1
        for i in range(dimx):
            for j in range(dimy):
                for k in range(dimz):
                    field_rnk = findRank(i, j, k)
                    x0, y0, z0 = x_y_z_list[field_rnk]
                    if tuples:
                        fld0[i, j, k] = np.array(fld_lst[field_rnk])[i - x0, j - y0, k - z0]
                    else:
                        fld0[i, j, k] = fld_lst[field_rnk][i - x0, j - y0, k - z0]
        return fld0
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
        data['mblocks_xyz'] = x_y_z_list
        data['mblocks_sxyz'] = sx_sy_sz_list
        # saving fields
        for f in range(nflds):
            fld_list = []
            for rnk in range(mpi_size):
                sx, sy, sz = sx_sy_sz_list[rnk]
                fld = np.array(struct.unpack("f" * sx * sy * sz, fileContent[read_ptr : read_ptr + 4 * sx * sy * sz]))
                read_ptr += 4 * sx * sy * sz
                fld = fld.reshape(sx, sy, sz) # then `fld[xi,yi,zi]` is the field at `xi,yi,zi`
                fld_list.append(fld)
            fld_list = globalizeFld(fld_list)
            data[variables[f]] = np.array(fld_list)
        # saving grid as a field
        fld_list = []
        for rnk in range(mpi_size):
            sx, sy, sz = sx_sy_sz_list[rnk]
            x0, y0, z0 = x_y_z_list[rnk]
            xyz_grid = [[[(i + x0, j + y0, k + z0) for k in range(sz)] for j in range(sy)] for i in range(sx)]
            fld_list.append(xyz_grid)
        fld_list = globalizeFld(fld_list, True)
        data['xyz'] = np.array(fld_list)
        if (len(fileContent) != read_ptr):
            print ("WRONG reading!")
    return data
# 2D usage example:
# ```
#   data = isolde.getFields("flds.tot.00000")
#   x_ = data['xyz'][:,:,0,0]
#   y_ = data['xyz'][:,:,0,1]
#   ex_ = data['ex'][:,:,0]
#   plt.pcolormesh(x_, y_, ex_) # <- 2D plot
# ```
