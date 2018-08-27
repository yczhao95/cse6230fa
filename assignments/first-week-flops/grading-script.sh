#!/bin/bash

[[ -z "$CSE6230_DIR" ]] && echo "CSE6230_DIR not defined, aborting" && exit 1
[[ -z "$GTUSER" ]] && echo "GTUSER not defined, aborting" && exit 1

module use ${CSE6230_DIR}/modulefiles
module unload cse6230
module load cse6230
module list

cd ${CSE6230_DIR}/assignments/first-week-flops/

[[ ! -z `git diff-files` ]] && echo "WARNING: unstaged changes!"

jupyter nbconvert --to=script first-week-flops.ipynb

# get the last qsub command in the file
last_qsub=`cat first-week-flops.sh | grep -n qsub | cut -f1 -d: | sort -nr | head -n1`
(head -${last_qsub} > first-week-flops-head.sh; cat > first-week-flops-compute.sh) < first-week-flops.sh

# Remove comments
sed -i -e '/^\s*qsub\b/ s/#.*//' first-week-flops-head.sh

# Make the script that gets submitted the compute node script
chmod u+x first-week-flops-head.sh
sed -i '1i #!/bin/bash\nset -x' first-week-flops-head.sh
sed -i -e '/^\s*qsub\b/ s/$/-x ${CSE6230_DIR}\/assignments\/first-week-flops\/first-week-flops-compute.sh/' first-week-flops-head.sh

# Make the compute node script executable
chmod u+x first-week-flops-compute.sh
sed -i '1i #!/bin/bash\ncd ${CSE6230_DIR}/assignments/first-week-flops\nmake --silent clean' first-week-flops-compute.sh

./first-week-flops-head.sh &> ./first-week-flops.${GTUSER}.grade.log

rm first-week-flops-head.sh
rm first-week-flops-compute.sh
