"""
author: Ian Beddows
email: ian.beddows@hhu.de
date: 10 Nov 2015

Author: Ian Beddows
Date: Nov 2015
Email: ian.beddows@hhu.de
"""
import numpy
import dadi
#~ from dadi import *
#=======================================================================
#-----------------------------------------------------------------------
#_________________________________F1____________________________________
#-----------------------------------------------------------------------
#=======================================================================

#split with symmetric migration between the two derived populations
def split_with_mig(params,(n1,n2), pts):
	TF,nu1F,nu2F,m12,m21=params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi , xx, TF, nu1F , nu2F, m12=m12, m21=m21)
	fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
	return fs


