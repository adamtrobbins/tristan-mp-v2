from tristanVis.alg import integrateFieldline2D
from tristanVis.tristanVis import Simulation
import tristanVis.aux as aux

class PulsarSimulation(Simulation):
  def loadData(self):
    import h5py
    import xarray as xr
    import numpy as np
    np.seterr(divide='ignore', invalid='ignore')
    with h5py.File(self._root + 'flds.tot.%05d' % self._step, 'r') as fields:
      with h5py.File(self._root + 'params.%05d' % self._step, 'r') as params:
        axes = ('z', 'y', 'x')
        self.data = xr.Dataset()
        self.data.attrs['t'] = params['timestep'][:][0]
        self.data.attrs['sx'], self.data.attrs['sy'], self.data.attrs['sz'] = (params['grd:mx0'][:][0], params['grd:my0'][:][0], params['grd:mz0'][:][0])
        xc = self.data.attrs['sx'] / 2
        yc = self.data.attrs['sy'] / 2
        zc = self.data.attrs['sz'] / 2
        self.data.attrs['CC'] = params['alg:c'][:][0]
        self.data.attrs['PERIOD'] = params['prb:psr_period'][:][0]
        self.data.attrs['RADIUS'] = params['prb:psr_radius'][:][0]
        self.data.attrs['OMEGA'] = 2 * np.pi / self.data.attrs['PERIOD']
        self.data.attrs['RLC'] = self.data.attrs['CC'] / self.data.attrs['OMEGA']
        self.data.attrs['sigma0'] = params['pls:sigma'][:][0]
        self.data.attrs['ppc0'] = params['pls:ppc0'][:][0]
        self.data.attrs['c_omp0'] = params['pls:c_omp'][:][0]
        self.data.attrs['qe'] = self.data.attrs['CC']**2 / (self.data.attrs['ppc0'] * self.data.attrs['c_omp0']**2)
        self.data.attrs['B0'] = self.data.attrs['CC']**2 * np.sqrt(self.data.attrs['sigma0']) / self.data.attrs['c_omp0']
        self.data.attrs['me'] = self.data.attrs['qe']
        self.data.attrs['nGJ'] = 2 * self.data.attrs['OMEGA'] * self.data.attrs['B0'] / (self.data.attrs['CC'] * self.data.attrs['qe'])

        b_sqr = (fields['bx'][:] * fields['bx'][:] + fields['by'][:] * fields['by'][:] + fields['bz'][:] * fields['bz'][:])
        self.data['jx'] = (axes, fields['jx'][:])
        self.data['jy'] = (axes, fields['jy'][:])
        self.data['jz'] = (axes, fields['jz'][:])

        self.data['b'] = (axes, np.sqrt(b_sqr))
        self.data['rho+'] = (axes, fields['dens2'][:])
        self.data['rho-'] = (axes, fields['dens1'][:])
        dens_tot = (self.data['rho+'] + self.data['rho-'])
        self.data['enrg+'] = (axes, fields['enrg2'][:])
        self.data['enrg-'] = (axes, fields['enrg1'][:])
        self.data['sigma'] = self.data.attrs['sigma0'] * (b_sqr / self.data.attrs['B0']**2) * (self.data.attrs['ppc0'] / (self.data['enrg+'] + self.data['enrg-']))
        self.data['gmean'] = (self.data['enrg+'] + self.data['enrg-']) / dens_tot
        self.data['de'] = self.data.attrs['c_omp0'] * np.sqrt(self.data['gmean']) * (self.data.attrs['ppc0'] / dens_tot)
        self.data['rL'] = self.data['gmean'] * (self.data.attrs['c_omp0'] / np.sqrt(self.data.attrs['sigma0'])) * (self.data.attrs['B0'] / np.sqrt(b_sqr))
        temp = (fields['bx'][:] * fields['ex'][:] + fields['by'][:] * fields['ey'][:] + fields['bz'][:] * fields['ez'][:]) / b_sqr
        self.data['e.b'] = (axes, np.nan_to_num(temp))
        r_cyl = np.sqrt((fields['xx'][:] - xc)**2 + (fields['yy'][:] - yc)**2)
        self.data['Omega-y'] = (axes, np.abs(((-fields['bz'][:] * fields['ex'][:] + fields['bx'][:] * fields['ez'][:]) / b_sqr) * self.data.attrs['RLC'] / r_cyl))
        self.data['rhoGJ'] = (axes, -2 * self.data.attrs['OMEGA'] * fields['bz'][:] / (self.data.attrs['CC'] * self.data.attrs['qe']))
        self.data['bx'] = ((axes), fields['bx'][:])
        self.data['by'] = ((axes), fields['by'][:])
        self.data['bz'] = ((axes), fields['bz'][:])

        self.data.coords['x'] = (('x'), (fields['xx'][:][0,0,:] - self.data.attrs['sx'] / 2) / self.data.attrs['RLC'])
        self.data.coords['y'] = (('y'), (fields['yy'][:][0,:,0] - self.data.attrs['sy'] / 2) / self.data.attrs['RLC'])
        self.data.coords['z'] = (('z'), (fields['zz'][:][:,0,0] - self.data.attrs['sz'] / 2) / self.data.attrs['RLC'])

  def drawData(self, savefig=None, fontsize=None, figsize=None):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    import numpy as np
    from matplotlib import rc

    if figsize is not None:
      dims = figsize
    else:
      if savefig is not None:
        dims = (40, 24)
      else:
        dims = (46, 24)
    if fontsize is not None:
      fontsize = fontsize
    else:
      if savefig is not None:
        fontsize = 25
      else:
        fontsize = 12

    super().drawData()
    rc('font',**{'size':fontsize})
    np.seterr(divide='ignore', invalid='ignore')
    # integrating fieldlines
    def stopIf(point):
      return (np.linalg.norm(point) < self.data.attrs['RADIUS'] / self.data.attrs['RLC']) or (np.linalg.norm(point) > (0.5 * (self.data.attrs['sx'] - 40) / self.data.attrs['RLC']))

    fieldlines_xz = []
    RR = self.data.attrs['RADIUS'] / self.data.attrs['RLC']
    thetas = np.linspace(0, 2*np.pi, 100)[:-1]
    xs = RR * np.cos(thetas)
    zs = RR * np.sin(thetas)
    for x, z in zip(xs, zs):
      fieldline = integrateFieldline2D([x, z], self.data.coords['x'].values, self.data.coords['z'].values, self.data['bx'].sel(y=0).values, self.data['bz'].sel(y=0).values, +1, stop_condition=stopIf)
      fieldlines_xz.append(fieldline)
      fieldline = integrateFieldline2D([x, z], self.data.coords['x'].values, self.data.coords['z'].values, self.data['bx'].sel(y=0).values, self.data['bz'].sel(y=0).values, -1, stop_condition=stopIf)
      fieldlines_xz.append(fieldline)

    fig = plt.figure(figsize=dims)

    nx = 4
    ny = 3
    nn = 1

    ax = plt.subplot(ny, nx, nn)
    im = self.data['rho-'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1, vmax=1e4), cmap='turbo', interpolation='gaussian')
    ax.set_aspect(1)
    fig.get_axes()[-1].axhline(self.data.attrs['nGJ'], lw=2.5, c='white')

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['rho-'].sel(z=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1, vmax=1e4), cmap='turbo', interpolation='gaussian')
    ax.set_aspect(1)
    fig.get_axes()[-1].axhline(self.data.attrs['nGJ'], lw=2.5, c='white')

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['rho+'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1, vmax=1e4), cmap='turbo', interpolation='gaussian')
    ax.set_aspect(1)
    fig.get_axes()[-1].axhline(self.data.attrs['nGJ'], lw=2.5, c='white')

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['rho+'].sel(z=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1, vmax=1e4), cmap='turbo', interpolation='gaussian')
    ax.set_aspect(1)
    fig.get_axes()[-1].axhline(self.data.attrs['nGJ'], lw=2.5, c='white')

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    (self.data['b'] / self.data.attrs['B0']).sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1e-2, vmax=1), cmap='jet', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['sigma'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=0.1, vmax=1e3), cmap='idl', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['gmean'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=1, vmax=1e4), cmap='fire', interpolation='gaussian')
    ax.set_aspect(1)

    for fieldline in fieldlines_xz:
      ax.plot(*fieldline.T, c='white', lw=0.5)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['Omega-y'].sel(y=0).plot.imshow(norm=mpl.colors.Normalize(vmin=0.5, vmax=1.5), cmap='jet', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['e.b'].sel(y=0).plot.imshow(norm=mpl.colors.Normalize(vmin=-0.1, vmax=0.1), cmap='bipolar', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['de'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=0.1, vmax=10), cmap='inferno', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    self.data['rL'].sel(y=0).plot.imshow(norm=mpl.colors.LogNorm(vmin=0.1, vmax=100), cmap='inferno', interpolation='gaussian')
    ax.set_aspect(1)

    nn += 1
    ax = plt.subplot(ny, nx, nn, sharex=ax, sharey=ax)
    im = ((self.data['rho+'] - self.data['rho-']) / self.data['rhoGJ']).sel(y=0).plot.imshow(norm=mpl.colors.SymLogNorm(vmin=-10, vmax=10, linthresh=0.1), cmap='bipolar', interpolation='gaussian')
    ax.set_aspect(1)
    im.colorbar.set_label('rho / rhoGJ')

    for ax in fig.get_axes()[0:-1:2]:
      ax.add_patch(Circle((0, 0), self.data.attrs['RADIUS'] / self.data.attrs['RLC'], zorder=200, edgecolor='magenta', facecolor='none', lw=2))

    plt.suptitle('t = ' + str(self.data.attrs['t'] / self.data.attrs['PERIOD']) + ' rotations', fontsize=int(fontsize*2));
    plt.tight_layout(rect=[0.04, 0.04, 0.94, 0.94])
    fig.subplots_adjust(hspace=0.2, wspace=0.3)
    self.saveFig(savefig)
