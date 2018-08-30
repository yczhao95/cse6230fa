
# coding: utf-8

# # Assignment 1: Performance Metrics & First Week Flop/s
# 
# This assignment has some questions that you need to answer with text, and some code that you need to write.
# 
# You should put all of you textual answers in this notebook: `Insert->Insert Cell Below` to create a new cell below
# the question, and `Cell->Cell Type->Markdown` to make it a cell for entering text.
# 
# You will test your code on the compute nodes of pace-ice, and that it also where we will evaluate it.
# Please complete the text portions when you are logged into a head node working locally, and leave the compute nodes for when you actually need them.

# ## Performance Metrics
# 
# In class we talked about the _strong-scaling efficiency_ of a parallel algorithm / machine pair: $H_f(P) = T_f(1) / (P T_f(P))$.
# 
# We then talked about the _weak-scaling efficiency_ of algorithm $f$ that can be applied to different problem sizes $N$: $E_f(N,P) = T_f(N/P,1) / T_f(N,P)$.
# 
# The question came up of how they are related to each other.
# 
# First, the notion of strong scaling doesn't have a concept of problem size, so let's add it: let's define
# 
# $$H_f(N,P) = T_f(N,1) / (P T_f(N,P)).$$
# 
# This is simply strong-scaling efficiency for each problem instance individually.
# 
# **Question 1 (1 pt):** Show that the relative order of strong and weak scaling efficiency can be related to the efficiency of the serial algorithm, that is, whether $T_f(N,1)$ as a function of $N$ exhibits superlinear or sublinear behavior.

# ## PACE-ICE
# 
# **Head node exercise 1 (1 pt):** What command should you run from a head node to see a list of all the compute nodes and their availability?

# In[1]:


get_ipython().run_cell_magic('bash', '', 'pace-check-queue coc-ice')


# Try it out: open up this notebook on a head node and compare the list you get to the [orientation slides](http://pace.gatech.edu/sites/default/files/pace-ice_orientation_0.pdf).  You'll see that it has grown, and they haven't updated the orientation slides.  We'll just have to find out what all these new nodes are for ourselves.

# ---
# ### A word on running jupyter on pace-ice:
# 
# As you've probably noticed, the screen refresh can be a bit laggy.  I'm looking for a solution for that.  In the mean time, know that you don't have to work directly in the notebook: you can work on you answers in the terminal, and then paste them into the notebook, as long as you're confident that they are correct.
# 
# ---

# **Head node exercise 2 (1 pt):** For the next questions, I need you to log in to compute nodes to find out about them, but you need to be able to specify which type of compute nodes you are accessing.
# 
# For each of the types of nodes that you see in the list of resources you've created, give me the `qsub` command to start an interactive job on that type of node, with the following requirements:
# 
# * The job should give you exclusive access to one node and all its cores and devices.
# * The job should let you pop open an X window (like a notebook) if you want to.
# * The job should begin in the CSE6230 directory.
# * The job should end after 30 minutes.

# In[13]:


# put one qsub command in this cell, and duplicate the cell for the others
qsub -I -X -q coc-ice -l nodes=1:ppn=4:gpus=1,walltime=00:30:00 -d $CSE6230_DIR


# ## What have we got to work with?
# 
# Now, we need to switch from a notebook running on the head node to one running on a compue node, so `File->Save and Checkpoint` this notebook and `File->Close and Halt` it.  (Now would also be a good time to `git add` and `git commit` changes to this file.)  Use one of your ineractive job scripts to connect to a compute node and run the notebook there.
# See you on the other side!
# 
# ---
# 
# Okay, you're running on the compute node.
# 
# **Compute node exercise 1 (1 pt):** Using bash scripting (`awk`, `grep`, `sed`) or any other tool you like (you could, e.g., write a python script in a separate file and call it, as long as you `git add` it), set the following variables so that the printout that follows is correct.  You script should be correct on any type of compute node.

# In[9]:


CPU_NAME=cat /proc/cpuinfo | grep -oP 'model name\s:\s\K.*' | head -1
CORE_COUNT=cat /proc/cpuinfo | grep 'model name' | wc -l
GPU_NAME=nvidia-smi --query-gpu=gpu_name --format=csv | grep -v name| uniq
GPU_COUNT=nvidia-smi --query-gpu=gpu_name --format=csv | grep -v name| wc -l


# In[7]:


echo "This nodes has ${CORE_COUNT} cores: its architecture is (Manufacturer, Product Id) ${CPU_NAME}"
if [[ ! $GPU_COUNT || $GPU_COUNT == 0 ]] ;  then
    echo "This node has no GPUs"
else
    echo "This node has ${GPU_COUNT} GPUs: its/their architecture is (Manufacturer, Product Id) ${GPU_NAME}"
fi


# **Compute node exercise 2 (1 pt):** After you have logged out of the compute node, use whatever resources published on the web you can find to estimate the peak single precision flop/s of this node (you only need to do this step for one of the types of nodes, not all of them).

# ## Flop/s fever
# 
# We've got to scratch that itch: we just want to go fast.  Okay, let's get it out of our system, and we'll look at more practical computations in future assignments.
# 
# You should choose one of the node types for this task.  Because this is more complex if multiple devices are involved
# **1 bonus point** is earned for choosing a node with GPUs.
# 
# **Compute node exercise 3 (2 pts):** The command below will compile and runs essentially the following computation:
# 
# ```C
# for (i = 0; i < N; i++) {
#   for (j = 0; j < T; j++) {
#     a[i] = a[i] * b + c;
#   }
# }
# ```
# And it will report the flop/s for the whole calculation.
# 
# Except that the array `a` will be spread out: `Nh` entries will be on the host and `Nd` entries will be on each of the devices.  Try to find values of `Nh`, `Nd`, and `T`, and (optionally) compiler optimization flags that give you the highest flop/s.  Things to consider:
# 
# - Try to make your whole computation run for about a second.
# - The time reported is the maximum time for any device: if one sits idle while the other finishes, it will rob you of flop/s.
# - I suggest looking at one type of device at a time: set one of `Nh` or `Nd` to zero.  Once you've found your best flop/s for that device, optimize the other, and then try to strike a balance.
# - Experiment with the merits of putting more weight on `Nh` and `Nd` vs more weight on `T`.
# - You can also choose to pass the option `Bs=X` to control the thread block size for the GPU, where `X` is a power of 2 between 64 and 2048.

# In[3]:


make run_fma_prof Nh=256 Nd=256 T=256 COPTFLAGS='-O -xHost' CUOPTFLAGS='-O' # modify this for peak flop/s


# **Compute Node Exercise 4 (2 pts):** Now let's see if we can make any transformations to the code to make a difference.
# 
# We will run the same program, but with fused multiply add loops that you have tried to optimize.  You should edit the files
# `fma_loop_host_opt.cu` and/or `fma_loop_dev_opt.c`: they start out exactly the same as the reference implementations used above.

# In[ ]:


diff fma_loop_host.c fma_loop_host_opt.c


# In[ ]:


diff fma_loop_dev.cu fma_loop_dev_opt.cu


# See if you can exploit vectorization, instruction level parallelism, and/or loop transformations to get a boost.

# In[ ]:


make run_fma_prof_opt Nh=256 Nd=256 T=256 COPTFLAGS='-O -xHost' CUOPTFLAGS='-O' # modify this for peak flop/s


# ## Submitting this work
# 
# **Workstation exercise 1 (1 pt):** When you have completed the rest of this assignment, `git add` the changes to this file, the source files you modified, and any scripts you added, and `git commit` them.
# 
# Assignments in this class are submitted by [transfering a copy][1] or your repository to the TA.
# 
# You do this by:
# 
# 1. Creating a duplicate repository on github with a very specific name.  The name should be `cse6230-YYY-gtusername`, where `YYY` is a short name for the assignment.  In this case, `YYY` is `a1`.  Make it a private repository.
# 
# 2. Push the version of your work that you want to be graded to the `master` branch of that repository.
# 
# 3. Request a transfer of ownership to the organization `cse6230`, team `Admin`.
# 
# If that seems like a lot to remember, I've made a script that does it for you.  If you have `GTUSER` defined as your userid in your environment, you can run the following script from the top directory:
# 
# ```
# ./utils/submit_assignment.sh master a1
# ```
# 
# This means that you intend to submit your `master` branch for assignment `a1`.  If you did your development for this assignment on a different branch, say `assignment-1`, you would run
# 
# ```
# ./utils/submit_assignment.sh assignment-1 a1
# ```
# 
# [1]: https://help.github.com/articles/transferring-a-repository-owned-by-your-personal-account/
# 
