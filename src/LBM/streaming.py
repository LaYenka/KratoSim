# this function compute the streaming step

from numpy import *

def streaming(self,fout,fin,v):
    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(roll(fout[i,:,:], v[i,0], axis=0),
                           v[i,1], axis=1 )
    return fin                       
