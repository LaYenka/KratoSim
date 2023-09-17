#
# Material Properties
#

import numpy as np
import math

class material:

  def struct_isotropic(self):
    E = 100.0
    v = 0.3
    C = E/(1.0+v)/(1.0-2.0*v) * np.array([[1.0-v, v, 0.0], [v, 1.0-v, 0.0], [0.0, 0.0, 0.5-v]])
    return E,v,C
