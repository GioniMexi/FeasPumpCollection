/**
 * @file fp_c_interface.h
 * @brief C interface to FP2.0
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2012
 */

#ifndef FP_C_INTERFACE_H
#define FP_C_INTERFACE_H

#include <ilcplex/cplex.h>

#ifdef __cplusplus
extern "C" {
#endif

int callFP(CPXENVptr env, CPXLPptr lp, int useFP2, int logOutput, double timeLimit, int* foundSol, double* foundObjval);

#ifdef __cplusplus
}
#endif

#endif /* FP_C_INTERFACE_H */
