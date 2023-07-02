# Feasibility Pump Collection

At this moment the code support the Feasibility Pump versions [[1]](#1),[[2]](#2),[[3]](#3),[[4]](#4).

The author of the current version of the code is Domenico Salvagnin dominiqs@gmail.com.



Compilation instructions
------------------------

This project uses CMake. Minimal steps (use the correct paths on your machine for CPLEX and XPRESS):

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -DCPLEX_ROOT_DIR=/opt/ilog/cos129/cplex -DXPRESSDIR=/opt/fico/xpressmp87 ..
- make -j12


Code overview
-------------

The main FP code is in feaspump.cpp. Interfaces to the supported LP solvers are in cpxmodel.cpp and xprsmodel.cpp.
Transformers (objects responsible for rounding a solution) are in transformers.cpp. Different strategies for sorting variables before rounding are in rankers.cpp.


Usage
-------------
```
$ .fp2 prob_file --config (-c) config_file
```

ToDo
-------------

- [ ] SoPlex/SCIP Interface
- [ ] Multiple Reference Vectors in the Projection Step
- [ ] Choices of Reference Vectors
- [ ] New Objective Scaling
- [ ] mRENS in Stage 3


References
------------
<a id="1">[1]</a> 
M. Fischetti, F. Glover, and A. Lodi. The feasibility pump. Mathematical
Programming, 104(1):91–104, 2005.

<a id="2">[2]</a> 
L. Bertacco, M. Fischetti, and A. Lodi. A feasibility pump heuristic for
general mixed-integer problems. Discrete Optimization, 4(1):63–76, 2007.

<a id="3">[3]</a> 
T. Achterberg and T. Berthold. Improving the feasibility pump. Discrete
Optimization,     4(1):77–86, 2007.

<a id="4">[4]</a> 
M. Fischetti and D. Salvagnin. Feasibility pump 2.0. Mathematical Program-
ming Computation, 1(2-3):201–222, 2009
