export  OMP_NUM_THREADS=4
export OMP_PROC_BIND=spread
export OMP_SCHEDULE=1

for N_p in 1 2 4 8 16; do
	  this_L=`echo "$N_p 0.333 20." | awk '{ print ($3 * $1^$2); }'`
	    this_T=`echo "$N_p 25600" | awk '{ print ($2 / ($1 * $1)); }'`
		  make runcloud NP=$(( 256*$N_p )) L=$this_L NT=$this_T PERF="perf stat"
done
