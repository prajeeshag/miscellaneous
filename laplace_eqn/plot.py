#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from subprocess import call

im=np.loadtxt('output.dat')
cs=plt.contourf(im)
plt.colorbar(cs)
plt.show()

