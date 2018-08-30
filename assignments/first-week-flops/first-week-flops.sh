
pace-check-queue coc-ice

qsub -I -X -q coc-ice -l nodes=1:ppn=4:gpus=1,walltime=00:30:00 -d $CSE6230_DIR

CPU_NAME=`cat /proc/cpuinfo | grep -oP 'model name\s:\s\K.*' | head -1`
CORE_COUNT=`cat /proc/cpuinfo | grep 'model name' | wc -l`
GPU_NAME=`nvidia-smi --query-gpu=gpu_name --format=csv | grep -v name| uniq`
GPU_COUNT=`nvidia-smi --query-gpu=gpu_name --format=csv | grep -v name| wc -l`

echo "This nodes has ${CORE_COUNT} cores: its architecture is (Manufacturer, Product Id) ${CPU_NAME}"
if [[ ! $GPU_COUNT || $GPU_COUNT == 0 ]] ;  then
    echo "This node has no GPUs"
else
    echo "This node has ${GPU_COUNT} GPUs: its/their architecture is (Manufacturer, Product Id) ${GPU_NAME}"
fi

make run_fma_prof Nh=256 Nd=256 T=256 COPTFLAGS='-O -xHost' CUOPTFLAGS='-O' # modify this for peak flop/s

diff fma_loop_host.c fma_loop_host_opt.c

diff fma_loop_dev.cu fma_loop_dev_opt.cu

make run_fma_prof_opt Nh=256 Nd=256 T=256 COPTFLAGS='-O -xHost' CUOPTFLAGS='-O' # modify this for peak flop/s
