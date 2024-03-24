# define lattice constants

from numpy import *

def D2Q9(self):
 ###### Lattice Constants #######################################################
 # Lattice velocity
 v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
 # Gauss Weights           
 t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])
 
 # directions
 col1 = array([0, 1, 2])
 col2 = array([3, 4, 5])
 col3 = array([6, 7, 8])
 
 return v,t,col1,col2,col3
