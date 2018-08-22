# CSE 6230 (Fall 2017, 3 Credit Hours) | High Performance Parallel Computing: Tools and Applications

## Contact

|          | **Instructor**                                      | **Teaching Assistant**                                      |
| ---      | :------------------                                 | :------------------                                         |
| *name*   | Toby Isaac                                          | Zhihang Liu                                                 |
| *email*  | [tisaac@cc.gatech.edu](mailto:tisaac@cc.gatech.edu) | [zliu354@gatech.edu](zliu354@gatech.edu)                    |
| *phone*  | [(404)385-5970](tel:+14043855970)                   |                                                             |
| *office* | [KACB 1302](https://goo.gl/38x6QB)                  | [KACB CSE Common Area](https://goo.gl/38x6QB)                                                           |
| *hours*  | Mondays 3-4, Fridays 10-11                          | Wednesdays 1-2

## Websites

We will use our [Canvas page][canvas] for announcements, discussion, quick questions, etc.

Your fellow students have a lot of expertise, and could benefit from yours.  Please take advantage of that with discussion on Canvas.

Most course materials (lecture notes and assignments) will be managed by the repository for this course:
[https://github.gatech.edu/cse6230fa18/cse6230](https://github.gatech.edu/cse6230fa18/cse6230).

**Note:** the GT GitHub server is only available from on campus or from behind
the campus [VPN].  This is a change from last year, and a bit of an annoyance.
However, our main computing resource this year, the [PACE Instructional
Cluster][PACE-ICE], is also behind the firewall, so having the course website
there is not a net annoyance.  If access to on campus resources from where you
would like to work is a problem, come talk to me.

[VPN]:https://oit.gatech.edu/services/end-point-computing/virtual-private-network-vpn
[PACE-ICE]:http://pace.gatech.edu/sites/default/files/pace-ice_orientation_0.pdf
[canvas]:https://gatech.instructure.com/courses/26393/

## Lectures

Tuesdays & Thursdays, 9:30--10:45 a.m., [College of Computing Building 101](https://goo.gl/maps/BDfZxqEkXw32)

## Description

So you've written an algorithm and coded it up.  Is your code any good?  Can it
be better?  What's the best it can be?

High Performance Computing (HPC) is about answering these
questions.  Answering these questions requires *performance analysis* of how
your algorithm matches up to the machine it is running on.  This course teaches
practical performance analysis to determine when and where code optimization
can help.

This course also teaches *how* to make your code better.  These days, making
your code better almost always means taking fullest advantage of the
*concurrency* in modern architectures.  This course is full of hands-on
experience with the software tools that turn efficient parallel algorithms into
efficient parallel programs.

## Objectives

I want you to be able to:

* *Design* high-performance code.  Good design requires good algorithms for
  the problem at hand -- the subject of other courses like [CSE 6220: Intro to
  HPC](http://cse6220.gatech.edu/sp17-oms/) -- but also good understanding of
  the computer system being used in the form of *performance models*.  This
  course will discuss practical performance modeling of current architectures,
  including multi-core, coprocessor (e.g. GPUs) and many-core designs.

* *Implement* high-performance code, drawing on a wide range of software
  engineering tools.  This course will have coding assignments that will
  familiarize you with what these tools can -- and can't -- do.

* *Analyze* high-performance code. This may require more sophisticated
  measurement than just the runtime on a given problem.  For performance to be
  reproducible and transferable, we must understand *why* it is or isn't fast.
  This course will present tools for measuring the behavior of different
  components of a computer system to identify bottlenecks and confirm whether
  the *implementation* matches the *design*.

## Prerequisite Courses and Skills

[CS 3120] is listed as a prerequisite on OSCAR.
Many of the concepts from that course listing (multi-threading, scheduling,
synchronization, and communication) will be relevant to this course, but
background references will be provided throughout.

This course can be thought of as the connection between [CSE 6220] and [CS
6290].  Neither is a prerequisite, but we will connect to some of the material
from each, and those who find themselves interested more in the algorithmic or
in the hardware aspects of this course can peruse the materials for those
courses for more depth.

Code examples and assignments from this course will be written in compiled
languages (C, C++, CUDA).  Most HPC computing systems run POSIX-y operating
systems, so basic comfort with command line tools is expected, and a
willingness to learn new tools will make life a lot easier.

**This is not a course for learning compiled languages for the first time:** In
this course, some assignments require you to write your code to match a given C
function signature like the following:

```C

int your_dgemm (int M, int N, int K,
                double *restrict A, double alpha,
                const double *restrict B,
                const double *restrict C);
```

If you come across problems in your code, such as:

* bad pointer dereferencing,
* running out of memory on the call stack,
* leaking memory on the heap,
* etc.,

and you do not feel confident in your ability to debug those problems for
yourself, you will have a difficult time in this class.  Ideally, most of your time
on assignments should be spent thinking about code *performance*, not code *correctness*.
An alternative course to consider would be [CSE 6010]: Computational Problem Solving.

[CS 3120]: https://oscar.gatech.edu/pls/bprod/bwckctlg.p_disp_course_detail?cat_term_in=201808&subj_code_in=CS&crse_numb_in=3210
[CSE 6220]: http://cse6220.gatech.edu/sp18-oms/
[CS 6290]: http://www.omscs.gatech.edu/cs-6290-high-performance-computer-architecture/
[CSE 6010]: https://oscar.gatech.edu/pls/bprod/bwckschd.p_disp_detail_sched?term_in=201808&crn_in=87078

## Course Material

The course will break down roughly into the following content modules.

1.

  * Introduction and course logistics
  * Machine model: Work/time, CPU & registers
  * Concepts:
    - scalability vs. efficiency
    - performance metrics
    - vectorization & instruction level parallelism
    - performance measurement
  * Tools:
    - compilers
    - hardware counters
    - basic CUDA

2.

  * Machine model: CPU(s) & (shared) memory
  * Concepts:
    - machine balance, arithmetic intensity & the roofline model
    - non-uniform memory access (NUMA) and pages
    - thread pinning
    - threads vs. processes
  * Tools:
    - STREAM benchmark
    - basic openMP
  * Applications:
    - N-body / molecular dynamics

3.

  * Machine model: CPUs, shared memory & cache
  * Concepts:
    - cache blocking & loop transformations
    - cache coherence
    - false sharing
  * Tools: 
    - pseudorandom number generators
    - openMP continued
    - CUDA continued
  * Applications:
    - stencils
    - shared-memory dense linear algebra

4.

  * Machine model: Nodes & networks
  * Concepts:
    - distributed memory models
    - collective operations
    - remote memory access
    - network topology
    - partitioning
  * Tools:
    - MPI
    - NCCL (? maybe)
    - HPCToolkit
  * Applications:
    - sorting
    - sparse linear algebra and graphs

5.

  * "Machine model": beyond the machine
  * Concepts:
    - Software engineering
    - Performance portability
    - Libraries & frameworks
  * Tools:
    - BLAS, LAPACK & other libraries
  * Additional topics as desired

## Assignments & Grading

|                     |     |
| ---                 | --- |
| Exercises (~1/week) | 1/3 |
| Projects (~3)       | 1/3 |
| Final Project       | 1/3 |

* **Exercises**: Exercises are mostly about practical experience.
  They are to be completed individually.
* **Projects**:  Projects will involve a combination of implementation and
  analysis to achieve high-performance.
  Students may work alone or in pairs.
* **Final Project**: The final project will involve each of the three main
  aspects of design, implementation, and analysis.
  Students may work alone or in pairs.


**Grades** will be earned according to the standard scale:

|     |          |
| --- | ---      |
| A   | 90--100% |
| B   | 80--90%  |
| C   | 70--80%  |
| D   | 60--70%  |
| F   | 0--60%   |

* **Passing** grades for Pass/Fail students are D or better.

* **Auditing** students are welcome to come to office hours to discuss
  exercises and projects; they may also submit code components of projects for
  automatic evaluation (e.g., we will put their submissions into the queue for
  benchmarking scripts).  They will not, however, receive graded feedback.

## Collaboration policy

For all assignments, students may collaborate through discussion, but all
coding and writing should be done by the project members.  Students who submit
unattributed material will be found in violation of the Honor code (see
Academic Integrity below).

## Late Policy and Resubmission

The first grade for any assignment will be based on the material available at
the due date: if nothing has been submitted, the first grade is 0.

An assignment may be resubmitted once after it is due to receive a second grade
of up to 85% of the points it would have received if it were submitted on time.

The final grade will be the maximum of the two grades.

Please be aware that grades for this course are due on December 17, 2018 at
12:00 PM.  To allow for sufficient time to grade all assignments, all
resubmissions that will be considered must be received by **December 13, 2018
at 11:59 PM, EST**.

Of course, extenuating circumstances and hardships will occur, which I will be
happy to discuss on a case-by-case basis.

## Textbooks and Reading

All required reading will be available online free to GT students.
Excellent resources include:

* [Professor Vuduc's lectures](http://vuduc.org/cse6230/)
* [Professor Chow's lectures](https://www.cc.gatech.edu/~echow/ipcc/hpc-course/)
* Georg Hager's lecture notes: follow the instructions under "Teaching
  material" on his [book's website](https://blogs.fau.de/hager/hpc-book)
* __Introduction to Parallel Computing__, by Grama et al.,
  Addison-Wesley, 2003. Available online via the Georgia Tech Library's
  website:
  [www](http://proquest.safaribooksonline.com/0-201-64865-2)
* Victor Eijkhout's [HPC book and course](https://bitbucket.org/VictorEijkhout/hpc-book-and-course/src/default/)

## Academic Integrity

Georgia Tech aims to cultivate a community based on trust, academic integrity,
and honor. Students are expected to act according to the highest ethical
standards.  For information on Georgia Tech's Academic Honor Code, please visit
[http://www.catalog.gatech.edu/policies/honor-code/](http://www.catalog.gatech.edu/policies/honor-code/)
or [http://www.catalog.gatech.edu/rules/18/](http://www.catalog.gatech.edu/rules/18/).

## Accomodations for Individuals with Disabilities

If you are a student with learning needs that require special accommodation,
contact the Office of Disability Services at [(404)894-2563](tel:+14048942563)
or
[http://disabilityservices.gatech.edu/](http://disabilityservices.gatech.edu/),
as soon as possible, to make an appointment to discuss your special needs and
to obtain an accommodations letter.  Please also e-mail me as soon as possible
in order to set up a time to discuss your learning needs.


