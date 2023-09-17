# save data 

import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import *

# Visualization of the velocity.
def plot_field(self,time,u):     
    if (time%1000==0):
        print('time...', time)
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        plt.savefig("vel.{0:04d}.Re.".format(time//100) + str(self.Re) + ".png")
