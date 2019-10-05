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

def split_no_mig(params,(n1,n2), pts):
	TF,nu1F,nu2F=params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi , xx, TF, nu1F , nu2F)
	fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
	return fs

def no_mig_size(params, ns, pts):
    """
    Split with no migration, then size change with no migration.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    """
    nu1a, nu2a, nu1b, nu2b, T1, T2 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)
    
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
	
#split with symmetric migration between the two derived populations
def split_with_mig(params,(n1,n2), pts):
	TF,nu1F,nu2F,m12,m21=params
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	phi = dadi.Integration.two_pops(phi , xx, TF, nu1F , nu2F, m12=m12, m21=m21)
	fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
	return fs

def asym_mig_size(params, ns, pts):
    """
    Split with different migration rates, then size change with different migration rates.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: Time in the past of split (in units of 2*Na generations)
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: Time of population size change.
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m12, m21=m21)
    
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

	
def sec_contact_asym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = dadi.Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
   
def sec_contact_asym_mig_size(params, ns, pts):
    """
    Split with no gene flow, followed by size change with asymmetrical gene flow.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between the secondary contact and present.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    """
    nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2 = params

    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)

    phi = dadi.Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m12, m21=m21)
    
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs	

