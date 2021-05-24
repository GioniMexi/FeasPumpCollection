# Feasibility Pump 3.0

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
