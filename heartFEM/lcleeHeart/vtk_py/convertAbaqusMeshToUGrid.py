########################################################################

import sys

from heartFEM.lcleeHeart.vtk_py.readAbaqusMesh import *
from heartFEM.lcleeHeart.vtk_py.writeUGrid     import *

########################################################################

name = sys.argv[1]

mesh = readAbaqusMesh(name + ".inp", "hex")
writeUGrid(mesh, name + ".vtk")
