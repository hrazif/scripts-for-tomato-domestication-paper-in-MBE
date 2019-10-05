# Numpy is the numerical library dadi is built upon
import numpy
from numpy import array

import dadi


# Estimate parameter uncertainties using the Godambe Information Matrix, to
# account for linkage in the data. 
import dadi.Godambe

# In demographic_models.py, we've defined a custom model for this problem
import demographic_models

# Load the data
data = dadi.Spectrum.from_file('/home/caicedo/Hamid/20190311_redoing_analyses_for_the_paper/7_dadi_dating/1_full_genome_SP_NECU_and_SLCSLL/same_sample_sizes_folded_with_masking/no_gene_flow[converged]/full_genome.fs')
data.mask[1,0] = True # mask singletons
data.mask[0,1] = True

#to mask invariable sites
data.mask[10,10] = True
data.mask[0,0] = True


ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [30,40,50]

all_boot = [dadi.Spectrum.from_file('./boots/boot_%s.fs' %ii) for ii in range(1,101)]
all_boot_fold = [_.fold() for _ in all_boot]



#comparing model 1  vs. model 2
#function, load from demographic_models.py
func = demographic_models.split_with_mig

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple best-fit params.  Since evalution is done with more
# complex model, need to insert zero value at corresponding
# parameter index in complex model. And we need to tell the LRT adjust function
# that the Xth parameter (counting from 0) is the nested one.

#for likelihood ratio testing
#TF,nu1F,nu2F,0,0
p_lrt = [0.08,0.07,0.14,0,0]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot_fold, p_lrt, data, nested_indices=[3,4], multinom=True)

ll_simple_model=-692152.7307
ll_complex_model=-364495.2053
D_adj = adj*2*(ll_complex_model - ll_simple_model)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))


# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))





#comparing model 1  vs. model 4
#function, load from demographic_models.py
func = demographic_models.no_mig_size

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

#for likelihood ratio testing
#complex model:nu1a, nu2a, nu1b, nu2b, T1, T2
#simple model: nu1,  nu2,  nu1,   nu2, T, 0

p_lrt = [0.07,0.14,0.07,0.14,0.08,0]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot_fold, p_lrt, data, nested_indices=[5], multinom=True)

ll_simple_model=-692152.7307
ll_complex_model=-627090.6832
D_adj = adj*2*(ll_complex_model - ll_simple_model)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))


# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting simpler model: {0:.4f}'.format(pval))





#comparing model 2  vs. model 3
#function, load from demographic_models.py
func = demographic_models.sec_contact_asym_mig

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


#for likelihood ratio testing
#complex model:nu1, nu2, m12, m21, T1, T2
#simple model: nu1, nu2,  m12, m21, 0, T

p_lrt = [0.16,0.29,0.72,0.96,0,2.53]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot_fold, p_lrt, data, nested_indices=[4], multinom=True)

ll_simple_model=-364495.2053
ll_complex_model=-184140.0616
D_adj = adj*2*(ll_complex_model - ll_simple_model)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))


# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting simpler model: {0:.4f}'.format(pval))






#comparing model 3  vs. model 5
#function, load from demographic_models.py
func = demographic_models.sec_contact_asym_mig_size

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


#for likelihood ratio testing
#complex model:nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2
#simple model: nu1,  nu2,  nu1,  nu2, m12, m21, T1, T2

p_lrt = [1.15,1.97,1.15,1.97,0.25,0.58,4.96,0.24]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot_fold, p_lrt, data, multinom=True)

ll_simple_model=-184140.0616
ll_complex_model=-81393.0875
D_adj = adj*2*(ll_complex_model - ll_simple_model)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))


# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting simpler model: {0:.4f}'.format(pval))





#comparing model 4  vs. model 5
#function, load from demographic_models.py
func = demographic_models.sec_contact_asym_mig_size

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple best-fit params.  Since evalution is done with more
# complex model, need to insert zero value at corresponding
# parameter index in complex model. And we need to tell the LRT adjust function
# that the Xth parameter (counting from 0) is the nested one.

#for likelihood ratio testing
#complex model:nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2
#simple model: nu1a, nu2a, nu1b, nu2b, 0, 0, T1, T2

p_lrt = [0.01,0.02,4.39,2.95, 0,0, 0.01,0.22]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot_fold, p_lrt, data, nested_indices=[4,5],multinom=True)


ll_simple_model=-627090.6832
ll_complex_model=-81393.0875
D_adj = adj*2*(ll_complex_model - ll_simple_model)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))


# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting simpler model: {0:.4f}'.format(pval))




#confidence interval of tau based on model 5

func = demographic_models.sec_contact_asym_mig_size

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

popt=[0.06,0.10,0.15,4.97,0.81,0.34,0.38,0.54]

uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot_fold, popt, data, multinom=True, return_GIM=True)

#to save GIM
numpy.savetxt('asym_mig_size_GIM_matrix.txt', uncerts[1], delimiter='\t')

#to get just the uncertainties for each parameter; the last number is for theta
uncertainties=numpy.sqrt(numpy.diag(numpy.linalg.inv(uncerts[1])))
numpy.savetxt('asym_mig_size_uncertainties.txt', uncertainties, newline='\t')

#to get the variance covariance matrix
var_cov=numpy.linalg.inv(uncerts[1]) 
numpy.savetxt('asym_mig_size_GIM_varcov.txt', var_cov, delimiter='\t')

# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.

#Nref: theta / (L * 4 * mu)
#tau: 2 * (TF) * Nref * g (generation)
theta=1181074.482
Nref=theta/(802027819*4*10**-8)
#time in years
TF=popt[6]+popt[7]
tau=2*TF*Nref

#sd of TF
#progation of uncertainties, see https://www.google.com/url?q=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FPropagation_of_uncertainty%23Example_formulas&sa=D&sntz=1&usg=AFQjCNEErbDzWUxwHaOrhWQSyNZqooLLIg

cov_T1_T2=var_cov[6,7]
TF_sd=numpy.sqrt((uncertainties[6]**2)+(uncertainties[7]**2)+2*var_cov[6,7])

#sd of Nref: sd of theta divided by (L * 4 * mu)
Nref_sd=uncertainties[8]/(802027819*4*10**-8)

#covariance of TF and Nref
# see http://www.kaspercpa.com/statisticalreview.htm

cov_T1_theta=var_cov[6,8]
cov_T2_theta=var_cov[7,8]
cov_TF_Nref=(cov_T1_theta + cov_T2_theta)/(802027819*4*10**-8)

#progation of uncertainties, see https://www.google.com/url?q=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FPropagation_of_uncertainty%23Example_formulas&sa=D&sntz=1&usg=AFQjCNEErbDzWUxwHaOrhWQSyNZqooLLIg
tau_sd=tau*numpy.sqrt((TF_sd/TF)**2+(Nref_sd/Nref)**2+(2*cov_TF_Nref/(TF*Nref)))

#95% confidence interval is +- 1.96*tau_sd
tau_sd




