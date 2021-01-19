import numpy as np

def DIVIDE(ax, ay, az, bx, by, bz):
  return (ax / bx, ay / by, az / bz)

def OVER(ax, ay, az, b):
  return (ax / b, ay / b, az / b)

def MULT(ax, ay, az, bx, by, bz):
  return (ax * bx, ay * by, az * bz)

def TIMES(ax, ay, az, b):
  return (ax * b, ay * b, az * b)

def CROSS(ax, ay, az, bx, by, bz):
  return (-(az * by) + ay * bz, az * bx - ax * bz, -(ay * bx) + ax * by)

def DOT(ax, ay, az, bx, by, bz):
  return ax * bx + ay * by + az * bz

def NORM(ax, ay, az):
  return np.sqrt(ax**2 + ay**2 + az**2)

def GAMMAtoBETA(gamma):
  return np.sqrt(1.0 - gamma**-2)

def BETAtoGAMMA(beta):
  return 1.0 / np.sqrt(1.0 - beta**2)

def UtoGAMMA(ux, uy, uz):
  return np.sqrt(1.0 + DOT(ux, uy, uz, ux, uy, uz))

def UtoBETA(ux, uy, uz):
  gamma = UtoGAMMA(ux, uy, uz)
  return (ux / gamma, uy / gamma, uz / gamma)

def VEC3(obj, name):
  return (obj[name + 'x'], obj[name + 'y'], obj[name + 'z'])

def boost_E(ex, ey, ez, bx, by, bz, ux, uy, uz):
  vx, vy, vz = UtoBETA(ux, uy, uz)
  gamma = UtoGAMMA(ux, uy, uz)
  beta_cross_B_x, beta_cross_B_y, beta_cross_B_z = CROSS(vx, vy, vz, bx, by, bz)
  beta_dot_e = DOT(vx, vy, vz, ex, ey, ez)
  dummy = gamma**2 / (gamma + 1.0)

  ex1 = gamma * (ex + beta_cross_B_x) - dummy * vx * beta_dot_e
  ey1 = gamma * (ey + beta_cross_B_y) - dummy * vy * beta_dot_e
  ez1 = gamma * (ez + beta_cross_B_z) - dummy * vz * beta_dot_e

  return (ex1, ey1, ez1)

def boost_B(bx, by, bz, ex, ey, ez, ux, uy, uz):
  vx, vy, vz = UtoBETA(ux, uy, uz)
  gamma = UtoGAMMA(ux, uy, uz)
  beta_cross_E_x, beta_cross_E_y, beta_cross_E_z = CROSS(vx, vy, vz, ex, ey, ez)
  beta_dot_b = DOT(vx, vy, vz, bx, by, bz)
  dummy = gamma**2 / (gamma + 1.0)

  bx1 = gamma * (bx - beta_cross_E_x) - dummy * vx * beta_dot_b
  by1 = gamma * (by - beta_cross_E_y) - dummy * vy * beta_dot_b
  bz1 = gamma * (bz - beta_cross_E_z) - dummy * vz * beta_dot_b

  return (bx1, by1, bz1)

def cart2sph(ax, ay, az, x, y, z, xc, yc, zc):
  rx = x - xc; ry = y - yc; rz = z - zc
  theta = np.arctan2(np.hypot(rx, ry), rz)
  phi = np.arctan2(ry, rx)
  ar = np.sin(theta) * (np.cos(phi) * ax + np.sin(phi) * ay) + np.cos(theta) * az
  atheta = np.cos(theta) * (np.cos(phi) * ax + np.sin(phi) * ay) - np.sin(theta) * az
  aphi = -np.sin(phi) * ax + np.cos(phi) * ay
  return (ar, atheta, aphi)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class VectorField:
  def __init__(self, fx, fy, fz):
    self._fx = fx
    self._fy = fy
    self._fz = fz
    self._sph = False
  def generateSpherical(self, xyz0):
    self._sph = True
    # ... meshgrid
    self._fr, self._ftheta, self._fphi = cart2sph(self._fx, self._fy, self._fz, *xyz, *xyz0)
  @property
  def x(self):
    return self._fx
  @property
  def y(self):
    return self._fy
  @property
  def z(self):
    return self._fz
  @property
  def r(self):
    assert self._sph
    return self._fr
  @property
  def theta(self):
    assert self._sph
    return self._ftheta
  @property
  def phi(self):
    assert self._sph
    return self._fphi

class Fields:
  __slots__ = ("dataset",)
  def __init__(self, ds):
    self.__class__.dataset.__set__(self, ds)
  def __getattr__(self, k):
    return getattr(self.dataset, k)
  def __setattr__(self, k, v):
    setattr(self.dataset, k, v)
  @property
  def e(self):
    return VectorField(getattr(self.dataset, 'ex'), getattr(self.dataset, 'ey'), getattr(self.dataset, 'ez'))
# class ScalarField:
#   def __init__(self, data, xyz):
#     self.data = data
#     self.xyz = xyz
#
# class VectorField:
#   def __init__(self, fx, fy, fz, xyz):
#     self._fx = ScalarField(fx, xyz)
#     self._fy = ScalarField(fy, xyz)
#     self._fz = ScalarField(fz, xyz)
#     self._sph = False
#   def generateSpherical(self, xyz0):
#     self._sph = True
#     self._fr, self._ftheta, self._fphi = cart2sph(self._fx, self._fy, self._fz, *xyz, *xyz0)
#   @property
#   def x(self):
#     return self._fx
#   @property
#   def y(self):
#     return self._fy
#   @property
#   def z(self):
#     return self._fz
#   @property
#   def r(self):
#     assert self._sph
#     return self._fr
#   @property
#   def theta(self):
#     assert self._sph
#     return self._ftheta
#   @property
#   def phi(self):
#     assert self._sph
#     return self._fphi
