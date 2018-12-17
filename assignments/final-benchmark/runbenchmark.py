import os
import string



#experiment on number of vertices(Nv = 2 ^ s)
for s in range(8,16):
	status = os.system(' ./graph500/seq-csr/seq-csr -s ' + str(s))

#experiment on edge factors
for e in range(8,16):
	status = os.system(' ./graph500/seq-csr/seq-csr -s 20 -e ' + str(e))
