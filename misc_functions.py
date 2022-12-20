"""Miscellaneous functions"""

import time

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
