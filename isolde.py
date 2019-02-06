import struct
import numpy as np

  # !--- PRTL.TOT.***** structure ----------------------------------!
  # ! HEADER:                                                      _
  # !   # of cpus.....................[4 bytes]                     |
  # !   # of species..................[4 bytes]                     |
  # !   # of variables................[4 bytes]                     |
  # !   variable names................[#var * 5 bytes]              |- disp_header
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
        mpi_size, nspec, nvars = struct.unpack("iii", fileContent[:4 * 3])
        variables = []
        variable_types = []
        read_ptr = 12
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
        data = {}
        for s in range(nspec):
            data[str(s)] = {}
            for i in range(nvars):
                (data[str(s)])[variables[i]] = np.array(struct.unpack(variable_types[i] * nprt[s], fileContent[read_ptr : read_ptr + nprt[s] * 4]))
                read_ptr += nprt[s] * 4
    return data
