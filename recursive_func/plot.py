#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from subprocess import call

while (True):
	x = int(raw_input('Enter the order of Legendre polynomial (Enter a negative number to exit) :'))

	if x < 0 :
		quit()

	Legend = []
	err = call(["./a.out", str(x)])

	if err != 0:
		print "Error while running a.out"
		quit()

	im=np.loadtxt('output.dat')
	plt.plot(im[:,0], im[:,1])
	Legend.append('P'+str(x))
	plt.legend(Legend)
	plt.show()

