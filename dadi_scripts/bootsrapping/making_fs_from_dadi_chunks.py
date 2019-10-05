#run this with python3
import dadi
from dadi import *
# scikit-learn bootstrap
from sklearn.utils import resample


#to make fs from dadi chunks
#number of bootstraps

n_chunks=85
for i in range(1,n_chunks):
	chunk = dadi.Misc.make_data_dict('../dadi_chunks/chunk%s.dadi' %i)
	pops = ['SP_NECU', 'SLC_SLL']
	pop_sizes = [10, 10]
	fs = dadi.Spectrum.from_data_dict(chunk, pops, pop_sizes, polarized=True) # fs object
	fs.mask[1,0] = True # mask singletons
	fs.mask[0,1] = True
	#to mask invariable sites
	fs.mask[10,10] = True
	fs.mask[0,0] = True
	#to save the fs file
	fs.to_file('fs_chunk%s.fs' %i) # save 'fs' bootstrap to file



#number of bootstraps
n_boot=100
for i in range(1,n_boot+1):
	chunk_nums = resample(range(1,n_chunks), replace=True, n_samples=n_chunks)
	fs = dadi.Spectrum.from_file('fs_chunk%s.fs' %chunk_nums[0])
	for j in range(1,n_chunks-1):
		fs2= dadi.Spectrum.from_file('fs_chunk%s.fs' %chunk_nums[j])
		fs=fs+fs2
	#to save the fs file
	fs.to_file('boot_%s.fs' %i) # save 'fs' bootstrap to file

	
