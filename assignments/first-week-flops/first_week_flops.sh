# This is the PBS script for running your first assignment
#PBS -N flops
#PBS -q coc-ice
#PBS -l walltime=00:05:00
#PBS -j oe
#PBS -o first_week_flops.out

cd $CSE6230_DIR/assignments/first-week-flops

# This is just a check to help you remember: if you have changes that aren't
# checked in, you may not be submitting all of your work!
git diff-files

# These variables are where you will put what you find about this nodes components
FWF_HOST_TYPE=
FWF_HOST_COUNT=
FWF_DEV_TYPE=
FWF_DEV_COUNT=

echo "This node's host architecture is [$FWF_HOST_TYPE], and has [$FWF_HOST_COUNT] cores"
if [ "x$FWF_DEV_COUNT" != "x" ]; then
  echo "This node's device architecture is [$FWF_DEV_TYPE], and has [$FWF_DEV_COUNT] devices"
else
  echo "This node does not have any accelerator devices"
fi

# You should change these values to try observe the highest joint flop rates
FWF_ND=256
FWF_NH=256
FWF_T=256

make clean
# This will run the baseline and optimized flop rate tests
# If you want to change any of the compiler flags that are used,
# you can do that by changing the command below to
#
# make run_fma_prof COPTFLAGS="-your -c -opt -flags" CUOPTFLAGS="-your -cu -opt -flags" Nd=$FWF_ND Nh=$FWF_NH T=$FWF_T
make run_fma_prof Nd=$FWF_ND Nh=$FWF_NH T=$FWF_T



