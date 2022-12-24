"""Miscellaneous functions"""

import time
import numpy as np

# timer for time benchmarking
time_benchmark = None

# Display time passed since last benchmark
def time_passed(init=False):
	global time_benchmark
	passed_time=0
	if time_benchmark!=None and not init:
		#passed_time = time.clock()-time_benchmark
		# No time.clock() since python 3.8. process_time should give CPU time. Otherise time.perf_counter() should be used, but on linux process_time() seems compatible with old time() call
		# passed_time = time.process_time()-time_benchmark
		passed_time = time.perf_counter()-time_benchmark
	#time_benchmark=time.clock()
	# time_benchmark=time.process_time()
	time_benchmark=time.perf_counter()
	return passed_time

def multiInterp2(x, xp, fp):
	"""1D interpolation over a 2D array
	from https://stackoverflow.com/questions/43772218/fastest-way-to-use-numpy-interp-on-a-2-d-array"""
	i = np.arange(x.shape[0])
	# i = np.arange(fp.shape[0])
	j = np.searchsorted(xp, x) - 1
	d = (x - xp[j]) / (xp[j + 1] - xp[j])
	print("aaa", i, j, fp.shape, d.shape)
	return (1 - d) * fp[i, j] + fp[i, j + 1] * d

# Check if there are no trees with the same name in the same files
def are_trees_unique(trees):
	names_files = []
	for tree in trees:
		names_files.append((tree.tree_name, tree.file_name))

	# Unique, if no repeating element in the names_files
	return len(set(names_files)) == len(names_files)
