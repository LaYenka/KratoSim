# -------------------------------------------------------------
# This class contains the shape function and 
# the associated gradients
# --------------------------------------------------------------

import numpy as np
import math
import sys 

class shapefun(object):

  def linear_2D(xi):
    """Shape functions for a 4-node, isoparametric element
		N_i(xi,eta) where i=[1,2,3,4]
		Input: 1x2,  Output: 1x4"""
    xi,eta = tuple(xi)
    N = [(1.0-xi)*(1.0-eta), (1.0+xi)*(1.0-eta), (1.0+xi)*(1.0+eta), (1.0-xi)*(1.0+eta)]
    """Gradient of the shape functions for a 4-node, isoparametric element.
		dN_i(xi,eta)/dxi and dN_i(xi,eta)/deta
		Input: 1x2,  Output: 2x4"""
    dN = [[-(1.0-eta),  (1.0-eta), (1.0+eta), -(1.0+eta)],
		  [-(1.0-xi), -(1.0+xi), (1.0+xi),  (1.0-xi)]]
    return 0.25 * np.array(N), 0.25 * np.array(dN)
  
  def linear_1D(xi):
    """Shape functions for a 2-node, isoparametric element
		N_i(xi) where i=[1,2]
		Input: 1x2,  Output: 1x2"""
    #xi = tuple(xi)
    N = [(1.0-xi), (1.0+xi)]
    """Gradient of the shape functions for a 2-node, isoparametric element.
		dN_i(xi)/dxi
		Input: 1x2,  Output: 1x2"""
    dN = [-1.0,  1.0]
    return 0.5*np.array(N), 0.5*np.array(dN)
  
  def beam_1D(xi):
    """Shape functions for a beam, isoparametric element
		N_i(xi) where i=[1,2]
		Input: 1x4,  Output: 1x4"""
    #xi = tuple(xi)
    N = [0.25*(xi**3 -3*xi +2), (xi**3-xi**2-xi+1),(xi**3-3*xi-2),(xi**3+xi**2-xi-1)]
    """Gradient of the shape functions for a 4-node, isoparametric element.
		dN_i(xi)/dxi
		Input: 1x4,  Output: 1x4"""
    dN = [0.25*(3*xi**2-3), (3*xi**2-2*xi-1), -0.25*(3*xi**2-3), (3*xi**2-2*xi-1)]
    # introduce the factor 2 to support gauss integration in the range [-1,1]
    # while formulas are written for the interval [-1,1]
    return 2*np.array(N), 2*np.array(dN)
  
	
