# measure.py

import primaldimer_py
import time

NUM_ITER = 1000

seq1 = {
    "CCAAACAAAGTTGGGTAAGGATAGATCAAT",
    "CCAAACAAAGTTGGGTAAGGATAGTTCAAT",
    "CCAACAAAGTTGGGTAAGGATAGATCAAT",
}
seq2 = {
    "ACTCCCATGGCATAGCTCCAAA",
    "ACTCCCATGGCATAGCTCCAGA",
    "CCTACTCCCATGGCATAACTCCAAA",
    "CTACTCCCATGGCATAGCTCCATA",
}
t0 = time.perf_counter()
for _ in range(NUM_ITER):
    primaldimer_py.do_pools_interact_py([*seq1], [*seq2], -26)
t1 = time.perf_counter()

took = (t1 - t0) / NUM_ITER
print(f"Took and avg of {took * 1000:.2f}ms per iteration")
