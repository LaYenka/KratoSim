# compute collision step

from numpy import *

def collision(self,omega,feq,fin,obstacle):
              # Collision step.
              fout = fin - omega * (fin - feq)

              # Bounce-back condition for obstacle.
              for i in range(9):
                    fout[i, obstacle] = fin[8-i, obstacle]
              return fout
