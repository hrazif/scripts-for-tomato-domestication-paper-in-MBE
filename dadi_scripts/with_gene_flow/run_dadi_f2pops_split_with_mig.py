#!/usr/bin/python
"""
author:	Ian Beddows		09/29/2015
contact	ian.beddows@hhu.de | beddowsi@msu.edu

Short description:

"""
import matplotlib
import random
matplotlib.use('Agg') # before import of pylab! use backend for graphing
import pylab
import csv
import dadi
from dadi import *
import numpy as np
import datetime
today = str(datetime.datetime.now())
import demographic_models
#=======================================================================
#fs file!
dadi_file = 'full_genome'
run_num = 1 # 
#=======================================================================
# variables
just_exit = 0 # just exit after 2d sfs
print2d_sfs = 1
print1d_sfs = 0
print2d_comp = 1
write2outfile = 1
verbose = 1
#
niter = 10000 # < niter iterations
n_indp = 10 # number of independent estimates of parameter values
n_bootsr = 100 # number of bootstrap estimates
pts_l = [50,60,70] # the grid point settings used for extrapolation.
#=======================================================================
#dd/fs in
pops = ['SP_NECU', 'SLC_SLL']


#these numbers are based on number of individuals per population multipled by 2, then by 0.9, to keep SNPs present in at least 90% of the accessions
pop_sizes = [18, 385] #



#@!dd = dadi.Misc.make_data_dict('%s.dadi' %dadi_file) # data dictionary
#@!fs = dadi.Spectrum.from_data_dict(dd, pops, pop_sizes, polarized=True) # fs object

#to check the number segregating sites
#@!fs.S()

#to save the fs file
#@!fs.to_file('%s.fs' %dadi_file) # save 'fs' object 2 file

#~ dd = dadi.Misc.make_data_dict('dadi_files/%s.dadi' %dadi_file) # data dictionary
#~ fs = dadi.Spectrum.from_data_dict(dd, pops, pop_sizes) # fs object
#~ fs.to_file('fs_files%s.fs' %dadi_file) # save 'fs' object 2 file
fs = dadi.Spectrum.from_file('%s.fs' %dadi_file)

fs.mask[1,0] = True # mask singletons
fs.mask[0,1] = True
#=======================================================================
"""
set the function and parameters:

"""
ns = fs.sample_sizes # sample_sizes is a method of the fs object
#~ fu = 'f2pops_mig'
#~ fu = 'f2pops_admix'
fu = 'split_with_mig'
#~ fu = 'f2pops_mig_bot2'
#=======================================================================
# print SFS
# 2d_sfs_fig
if print2d_sfs:
	tmp = pylab.figure()
	dadi.Plotting.plot_single_2d_sfs(fs,vmin=0.0001)
	tmp.savefig('2d_sfs_fig/%s-%s_%s_2d_sfs.pdf' %(pops[0],pops[1],dadi_file) )
if print1d_sfs:
	tmp = pylab.figure()
	dadi.Plotting.plot_1d_fs(fs)
	tmp.savefig('1d_fig/%s_%s_1d_fs.pdf' %(pops[0],dadi_file) )
if just_exit:
	import sys
	sys.exit()
	# exit!
#=======================================================================
'''
set parameters dependent on function
'''

if fu=='split_with_mig':
	func = demographic_models.split_with_mig
	params = ['TF','nu1F','nu2F', 'm12', 'm21']
	upper_bound = [5, 5, 5,10,10]
	lower_bound = [1e-2, 1e-2, 1e-2,1e-2,1e-2]#the lower bound should not be zero, because log of zero is NA.
	_p0 = [0.3161,0.2001,0.7687,0.0039,0.9784]#this is the highest likelihood estimates that I use as starting points for BS runs.


#=======================================================================
# define function for extrapolation
func_name = func.__name__
func_ex = dadi.Numerics.make_extrap_log_func(func)
#=======================================================================
# write to outfile
if write2outfile: 
	outfile = 'out/%s-%d_sfsFILE-%s.csv' %(func_name,run_num,dadi_file)
	f = open(outfile,'a') # filehandle
	w = csv.writer(f, lineterminator='\n') # writer object 'w'
	params.insert(0,'independent estimate number')
	params.insert(1,'maximum log composite likelihood')
	params.insert(2,'theta')
	f.write("\n\n# %s\n" %today)
	f.write("# %s\n" %dadi_file)
	f.write("# Function: %s\n" %func_name)
	f.write("# Pop: %s\n" %pops)
	f.write("# %i independent estimates\n" %n_indp)
	f.write("# < %i iterations\n" %niter)
	w.writerow(params)
#=======================================================================
'''
compute n_indp estimates of the parameter values (p0)


'''
for i in range(1,n_indp+1):
	p0=[]
	p0.append(random.uniform(1e-2,5))
	p0.append(random.uniform(1e-2,5))
	p0.append(random.uniform(1e-2,5))
	p0.append(random.uniform(1e-2,10))
	p0.append(random.uniform(1e-2,10))
	p0 = dadi.Misc.perturb_params(p0, fold=0, upper_bound=upper_bound, lower_bound=lower_bound)


	popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, # estimate
									   lower_bound=lower_bound,
									   upper_bound=upper_bound,
									   verbose=1, maxiter=niter)
	model = func_ex(popt, ns, pts_l) # Calculate the best-fit model AFS.
	ll_model = dadi.Inference.ll_multinom(model, fs) # Likelihood of the data given the model AFS.
	theta = dadi.Inference.optimal_sfs_scaling(model, fs)

	if print2d_comp:
		tmp = pylab.figure()
		dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=20,
										pop_ids =(pops))
										
		tmp.savefig('2d_comp_fig/%d_%s_%s_%s-%s.pdf' %(i,dadi_file,func_name,pops[0],pops[1]))

	if write2outfile:	# write to outfile
		out = np.array(popt).tolist()
		out.insert(0,i)
		out.insert(1,ll_model)
		out.insert(2,theta)
		out = ['%.4f' % i for i in out]
		w.writerow(out)

	if verbose:
		print('Best-fit parameters: {0}'.format(popt))
		print('Maximum log composite likelihood: {0}'.format(ll_model)) # The optimal value of theta given the model.
		print('Optimal value of theta: {0}'.format(theta))
#End

f.write('\n\n# %d bootstrap estimates\n' %n_bootsr)
#=======================================================================
# bootstrap
for i in range(1,n_bootsr+1):
	p0 = _p0
	p0 = dadi.Misc.perturb_params(p0, fold=0, upper_bound=upper_bound, lower_bound=lower_bound)

	fs_b = fs.sample()

	popt = dadi.Inference.optimize_log(p0, fs_b, func_ex, pts_l, # estimate
									   lower_bound=lower_bound,
									   upper_bound=upper_bound,
									   verbose=len(p0), maxiter=niter)
	model = func_ex(popt, ns, pts_l) # Calculate the best-fit model AFS.
	ll_model = dadi.Inference.ll_multinom(model, fs_b) # Likelihood of the data given the model AFS.
	theta = dadi.Inference.optimal_sfs_scaling(model, fs_b)
	
	if write2outfile:	# write to outfile
		out = np.array(popt).tolist()
		out.insert(0,i)
		out.insert(1,ll_model)
		out.insert(2,theta)
		out = ['%.4f' % i for i in out]
		w.writerow(out)
	
	if verbose:
		print('Best-fit parameters: {0}'.format(popt))
		print('Maximum log composite likelihood: {0}'.format(ll_model)) # The optimal value of theta given the model.
		print('Optimal value of theta: {0}'.format(theta))







# Plot a comparison of the resulting fs with the data.
#~ pylab.figure(1)
#~ dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=3,
                                    #~ pop_ids =('Schi','Sper'))
#~ pylab.savefig('prior_onegrow_mig.pdf')

#~ w.close() # close csv writer
f.close() # close filehandle
'''
------------------------------------------------------------------------
likelihood-ratio test comparing models 2 follow
------------------------------------------------------------------------
'''



