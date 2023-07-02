/**
 * @file pdlpmodel.cpp
 * @brief Implementation of MIPModelI for PDLP
 *
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
 */

#include "feaspump/pdlpmodel.h"
#include <signal.h>
#include <cstring>
#include <climits>
#include <numeric>
#include <fmt/format.h>
#include <cassert>
#include <unistd.h>
#include <utils/maths.h>
#include <utils/floats.h>
#include <utils/timer.h>

int PDLPModel::numGeneralInts()
{
	// col types are obtained from the cloned problem
	DOMINIQS_ASSERT(scip);
	SCIP_VAR ** vars;
	int nvars;
	int nints;
	SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
	nints = 0;

	for (int idx = 0; idx <= nvars; idx++)
	{
		if (SCIPvarIsIntegral(vars[idx]))
		{
			if (!SCIPvarIsBinary(vars[idx])) nints++;
		}
	}
	return nints;	
}

double PDLPModel::getPdlpTolerance()
{
	return params.mutable_termination_criteria()->eps_optimal_relative();
}

void PDLPModel::setPdlpTolerance(double tol, double decreaseFactor)
{
	tol = std::max(tol * decreaseFactor, 1e-8);
	params.mutable_termination_criteria()->set_eps_optimal_absolute(tol);
	params.mutable_termination_criteria()->set_eps_optimal_relative(tol);
	// params.mutable_termination_criteria()->set_eps_primal_infeasible(tol);
	// params.mutable_termination_criteria()->set_eps_dual_infeasible(tol);
}

SCIP_STAGE PDLPModel::getProbStage()
{
	DOMINIQS_ASSERT(scip);
	return SCIPgetStage(scip);
}

PDLPModel::PDLPModel()
{
	// initialize scip problem
	DOMINIQS_ASSERT(SCIPcreate(&scip));
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_PROBLEM));
	DOMINIQS_ASSERT(SCIPincludeDefaultPlugins(scip));
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
	// initialize PDLP
	solver = MPSolver::CreateSolver( "PDLP" ) ;
}


PDLPModel::PDLPModel(SCIP* _prob, bool _ownProb) : scip(_prob), ownProb(_ownProb)
{
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(solver);
}

/* Destroyer */
PDLPModel::~PDLPModel()
{
	// ToDo
}


/* Read/Write */
void PDLPModel::readModel(const std::string& filename)
{
	// read problem from file
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(SCIPreadProb(scip, filename.c_str(), NULL));
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_PROBLEM));
	DOMINIQS_ASSERT(SCIPtransformProb(scip));
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_TRANSFORMED));

}

void PDLPModel::writeModel(const std::string& filename, const std::string& format) const
{
    // ToDo
}

void PDLPModel::writeSol(const std::string& filename) const
{
    // ToDo
}


/* Solve */
double PDLPModel::lpopt(char method, bool decrease_tol, bool initial)
{
	// solve by PDLP and update the result_status
	DOMINIQS_ASSERT(solver);
	if (verbosity) solver->EnableOutput();
	else solver->SuppressOutput();
	absl::lts_20211102::StatusOr<operations_research::pdlp::QuadraticProgram> qp;
	double current_tol = getPdlpTolerance();
	setPdlpTolerance(current_tol, pdlpTolDecreaseFactor);
	if(decrease_tol)
	{
		setPdlpTolerance(current_tol, 0.1);
	}
	operations_research::MPModelProto output_model;

	dominiqs::StopWatch chrono;
	chrono.start();

	solver->ExportModelToProto(&output_model);
    qp = operations_research::pdlp::QpFromMpModelProto(output_model, true, false);
	DOMINIQS_ASSERT(qp.ok());
	operations_research::pdlp::QuadraticProgram real_qp = *qp;
	double time_model = chrono.getElapsed();
	if (pdlpWarmStart && !initial)
	{
		initial_solution = {output.primal_solution, output.dual_solution};
		initial_solution->primal_solution.resize(ncols());
		initial_solution->dual_solution.resize(nrows());
		// for the added cols and rows we initialize by zero since they may be different cols and rows
		// than those from the previous iteration
		for (int i = numRowsPresolved; i < nrows(); i++)
		{
			initial_solution->dual_solution[i] = 0; 
		}
		for (int i = numColsPresolved; i < ncols(); i++)
		{
			initial_solution->primal_solution[i] = 0; 
		}
		output = operations_research::pdlp::PrimalDualHybridGradient(
		real_qp, params, initial_solution);	
	}
	else output = operations_research::pdlp::PrimalDualHybridGradient(
         real_qp, params);

	result_status = output.solve_log.termination_reason();
	return time_model;
}

/* Solve mip */
void PDLPModel::mipopt()
{ 
    // ToDo
}

/* presolve mip */
void PDLPModel::presolve()
{
	// presolve the problem
	DOMINIQS_ASSERT(scip);
	
	dominiqs::StopWatch chrono;
	chrono.start();

	SCIPsetIntParam(scip, "display/verblevel", 0);

	SCIPsetBoolParam(scip,"constraints/linear/upgrade/indicator" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/knapsack" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/logicor" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/setppc" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/varbound" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/xor" , false);
	std::cout << "SCIP presolve" << std::endl;
	DOMINIQS_ASSERT(SCIPpresolve(scip));
	chrono.stop();
	std::cout << "SCIP presolve time: " << chrono.getElapsed() << std::endl;
	std::cout << "SCIP presolve done" << std::endl;

}

// Postsolve the problem
void PDLPModel::postsolve()
{
	// not needed since we can access original vars and constraints in scip
	if (!stageLP)
	{
		postsolved = true;
		DOMINIQS_ASSERT(scip);
	}
}

/* get solution vector in the original space */
std::vector<double> PDLPModel::postsolveSolution(const std::vector<double>& preX) const
{
	// create sol for the presolved MIP
	SCIP_SOL *sol;
	SCIPcreateSol(scip,  &sol, NULL);
	SCIP_VAR **vars;

	int nvars = SCIPgetNVars(scip);
	// SCIPallocBufferArray(scip, &vars, nvars);

	SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
	// set solution values for the variables
	for (int idx = 0; idx < nvars; idx++)
	{
		SCIPsetSolVal(scip, sol, vars[idx], preX[idx]);
	}
	// retransform sol to original space
	SCIPretransformSol(scip, sol);
	// get variables in original space
	SCIP_VAR **varsOrig;
	nvars = SCIPgetNOrigVars(scip);
	SCIPgetOrigVarsData(scip, &varsOrig, &nvars, NULL,NULL,NULL,NULL);
	// origX is the solution vector in the original space
	std::vector<double> origX(nvars, 0);
	for ( int i = 0; i < nvars; i++)
	{
		origX[i] = SCIPgetSolVal(scip, sol, varsOrig[i]);
	}
	SCIPfreeSol(scip, &sol);
	return origX;
}


/* Get solution */
double PDLPModel::objval() const
{
	// objective value from PDLP
	DOMINIQS_ASSERT(solver);
	std::vector<double> coefs(ncols(), 0);
	objcoefs(coefs.data(), 0, -1);
	double objval = 0.0;
	for (int i = 0; i < ncols(); i++) objval += (coefs[i] * output.primal_solution(i));
	return objval;
}


void PDLPModel::sol(double* x, int first, int last) const
{
	if (last == -1)  last = ncols()-1;
	if (stageLP)
	{
		for( int i = 0; i <= last; i++ )
	  	{

			x[i] = output.primal_solution(i);
	  	}
	}
}


bool PDLPModel::isPrimalFeas() const
{
	// ToDo check this once more
	// In case of PDLP we do not care if the solution is feasible or not
	// Therefore whenever the algorithm terminates without an error we return true
	if (stageLP)
	{
		DOMINIQS_ASSERT(solver);
		if ( result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_OPTIMAL || 
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_PRIMAL_INFEASIBLE ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_DUAL_INFEASIBLE ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_TIME_LIMIT ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_ITERATION_LIMIT ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_INTERRUPTED_BY_USER ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_NUMERICAL_ERROR ||
			 result_status == operations_research::pdlp::TerminationReason::TERMINATION_REASON_PRIMAL_OR_DUAL_INFEASIBLE 		 
			 )
		{
			return 1;
		}
		// Other cases:
		else
		{	
			return 0;
		}
    }
	// ToDo for stage 3...
	return 1;
}

/* Parameters */
void PDLPModel::handleCtrlC(bool flag)
{
	// ToDo
}


bool PDLPModel::aborted() const
{
	// ToDo
    return false;
}


void PDLPModel::seed(int seed)
{
	// ToDo
}


void PDLPModel::logging(bool log)
{
	// For SCIP
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(solver);
	verbosity = log;
	if (log)  SCIPsetIntParam(scip, "display/verblevel", 4);
	else      SCIPsetIntParam(scip, "display/verblevel", 0);
	// For PDLP
	if (log) solver->EnableOutput();
	else solver->SuppressOutput();

}

int PDLPModel::intParam(IntParam which) const
{
	int value;
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		switch(which)
		{
			case IntParam::Threads:
				SCIPlpiGetIntpar(lpi, SCIP_LPPAR_THREADS, &value);
				break;
			case IntParam::SolutionLimit:
				break;
			case IntParam::NodeLimit:
				break;
			case IntParam::IterLimit:
				SCIPlpiGetIntpar(lpi, SCIP_LPPAR_LPITLIM, &value);
				break;
			default:
				throw std::runtime_error("Unknown integer parameter");
		}

	}
	else
	{
		DOMINIQS_ASSERT(scip);
		switch(which)
		{
			case IntParam::Threads:
				SCIPgetIntParam(scip, "lp/threads", &value);
				break;
			case IntParam::SolutionLimit:
				SCIPgetIntParam(scip, "limits/solutions", &value);

				break;
			case IntParam::NodeLimit:
				SCIPgetIntParam(scip, "limits/totalnodes", &value);
				break;
			case IntParam::IterLimit:
				// not needed atm
				break;
			default:
				throw std::runtime_error("Unknown integer parameter");
		}
	}
	return value;
}


void PDLPModel::intParam(IntParam which, int value)
{
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		switch(which)
		{
			case IntParam::Threads:
				SCIPlpiSetIntpar(lpi, SCIP_LPPAR_THREADS, value);
				break;
			case IntParam::SolutionLimit:
				break;
			case IntParam::NodeLimit:
				break;
			case IntParam::IterLimit:
				params.mutable_termination_criteria()->set_iteration_limit(value);
				SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, value);
				break;
			case IntParam::PdlpWarmStart:
				if (value > 0) pdlpWarmStart = true;
				else pdlpWarmStart = false;
				break;
			default:
				throw std::runtime_error("Unknown integer parameter");
		}
	}
	else
	{
		DOMINIQS_ASSERT(scip);
		switch(which)
		{
			case IntParam::Threads:
				SCIPsetIntParam(scip, "lp/threads", value);
				break;
			case IntParam::SolutionLimit:
				SCIPsetIntParam(scip, "limits/solutions", value);

				break;
			case IntParam::NodeLimit:
				SCIPsetIntParam(scip, "limits/totalnodes", value);
				break;
			case IntParam::IterLimit:
				// not needed atm
				break;
			case IntParam::PdlpWarmStart:
				if (value > 0) pdlpWarmStart = true;
				break;
			default:
				throw std::runtime_error("Unknown integer parameter");
		}
	}

}


double PDLPModel::dblParam(DblParam which) const
{
	// ToDo for PDLP
	double value = 0.0;
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		switch(which)
		{
			case DblParam::TimeLimit:
				SCIPlpiGetRealpar(lpi, SCIP_LPPAR_LPTILIM, &value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPlpiGetRealpar(lpi, SCIP_LPPAR_FEASTOL, &value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPlpiGetRealpar(lpi, SCIP_LPPAR_FEASTOL, &value);
				break;
			case DblParam::PdlpTolerance:
				value = pdlpTol;
				break;
			case DblParam::PdlpToleranceDecreaseFactor:
				value = pdlpTolDecreaseFactor;
				break;
			default:
				throw std::runtime_error("Unknown double parameter");
		}
	}

	else
	{
		DOMINIQS_ASSERT(scip);
		switch(which)
		{
			case DblParam::TimeLimit:
				SCIPgetRealParam(scip, "limits/time", &value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPgetRealParam(scip, "numerics/feastol", &value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPgetRealParam(scip, "numerics/feastol", &value);
				break;
			case DblParam::PdlpTolerance:
				value = pdlpTol;
				break;
			case DblParam::PdlpToleranceDecreaseFactor:
				value = pdlpTolDecreaseFactor;
				break;
			default:
				throw std::runtime_error("Unknown double parameter");
		}
	}

	return value;
}


void PDLPModel::dblParam(DblParam which, double value)
{
	// ToDo for PDLP
	if (stageLP)
	{
		switch(which)
		{
			case DblParam::TimeLimit:
				// time limit for PDLP is set in milliseconds
				if (solver != NULL) solver->set_time_limit((int) value * 1000);
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LPTILIM, value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, value);
				break;
			case DblParam::PdlpTolerance:
				pdlpTol = value;
				setPdlpTolerance(pdlpTol, 1.0);
				break;
			case DblParam::PdlpToleranceDecreaseFactor:
				pdlpTolDecreaseFactor = value;
				break;
			default:
				throw std::runtime_error("Unknown double parameter");
		}
	}
	else
	{
		switch(which)
		{
			case DblParam::TimeLimit:
				// time limit for PDLP is set in milliseconds
				if (solver != NULL) solver->set_time_limit((int) value * 1000);
				SCIPsetRealParam(scip, "limits/time", value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPsetRealParam(scip, "numerics/feastol", value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPsetRealParam(scip, "numerics/feastol", value);
				break;
			case DblParam::PdlpTolerance:
				pdlpTol = value;
				setPdlpTolerance(pdlpTol, 1.0);
				break;
			case DblParam::PdlpToleranceDecreaseFactor:
				pdlpTolDecreaseFactor = value;
				break;
			default:
				throw std::runtime_error("Unknown double parameter");
		}

	}

}

int PDLPModel::intAttr(IntAttr which) const
{
	// ToDo
	int value = 0;
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		switch(which)
		{
			case IntAttr::Nodes:
				break;
			case IntAttr::NodesLeft:
				break;
			case IntAttr::BarrierIterations:
				break;
			case IntAttr::SimplexIterations:
				SCIPlpiGetIterations(lpi, &value);
			case IntAttr::PDLPIterations:
				DOMINIQS_ASSERT(solver);
				value = output.solve_log.iteration_count();
				break;
			default:
				throw std::runtime_error("Unknown integer attribute");
		}
	}
    return value;
}


double PDLPModel::dblAttr(DblAttr which) const
{
	// not needed
    return 0.0;
}

void PDLPModel::terminationReason(std::string& reason)
{
	// solve by PDLP and update the result_status
	DOMINIQS_ASSERT(solver);
	switch(result_status)
	{
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_OPTIMAL:
			reason = "OPTIMAL";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_TIME_LIMIT:
			reason = "TIME_LIMIT";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_ITERATION_LIMIT:
			reason = "ITERATION_LIMIT";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_UNSPECIFIED: 
			reason = "UNSPECIFIED";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_PRIMAL_INFEASIBLE:
			reason = "PRIMAL_INFEASIBLE";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_DUAL_INFEASIBLE:
			reason = "DUAL_INFEASIBLE";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT:
			reason = "KKT_MATRIX_PASS_LIMIT";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_INTERRUPTED_BY_USER:
			reason = "INTERRUPTED_BY_USER";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_NUMERICAL_ERROR:
			reason = "NUMERICAL_ERROR";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_INVALID_PROBLEM: 
			reason = "INVALID_PROBLEM";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_INVALID_PARAMETER:
			reason = "INVALID_PARAMETER";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_OTHER:
			reason = "OTHER";
			break;
		case operations_research::pdlp::TerminationReason::TERMINATION_REASON_PRIMAL_OR_DUAL_INFEASIBLE:
			reason = "PRIMAL_OR_DUAL_INFEASIBLE";
			break;
		default:
			reason = operations_research::pdlp::TerminationReason_Name(result_status);
			break;
	}
}


/* Access model data */
int PDLPModel::nrows() const
{
	// ToDo (check this when called in feaspump)
	DOMINIQS_ASSERT(solver);
    int ret;
	if (!isClone)
	{
		if (postsolved) ret = SCIPgetNOrigConss(scip);
		else ret = SCIPgetNConss(scip);
	}
	else
	{
		if (stageLP)
		{
			// DOMINIQS_ASSERT(SCIPlpiGetNRows(lpi, &ret));
			// ToDo (this is a walkround since nrows() is called at a point where solver is an empty problem)
			// Happens because we still do not have a way to just delete auxiliary cols and rows
			ret = solver->NumConstraints() != 0 ? solver->NumConstraints() : numRowsPresolved;
		}
		else
		{
			if (postsolved) ret = nconsStartStage3;
			else ret = SCIPgetNConss(scip);
		}
	}
	return ret;
}


int PDLPModel::ncols() const
{
	// ToDo (check this when called in feaspump)
    int ret;
	DOMINIQS_ASSERT(solver);
	if (!isClone)
	{
		if (postsolved) ret = SCIPgetNOrigVars(scip);
		else ret = SCIPgetNVars(scip);
	}
	else
	{
		if (stageLP)
		{
			// DOMINIQS_ASSERT(SCIPlpiGetNCols(lpi, &ret));
			// ToDo (this is a walkround since nrows() is called at a point where solver is an empty problem)
			// Happens because we still do not have a way to just delete auxiliary cols and rows
			ret = solver->NumVariables() != 0 ? solver->NumVariables() : numColsPresolved;

		}
		else
		{
			if (postsolved) ret = nvarsStartStage3;
			else ret = SCIPgetNVars(scip);
		}
	}
	return ret;
}


int PDLPModel::nnz() const
{
	// ToDo
    int ret;
	if (stageLP) DOMINIQS_ASSERT(SCIPlpiGetNNonz(lpi, &ret));
	else ret = SCIPgetNNZs(scip);
	return ret;
}

double PDLPModel::objOffset() const
{
	// ToDo
	double ret = 0.0;
	if (isClone) ret = SCIPgetOrigObjoffset(scip);
	else
	{
		if (postsolved) ret = SCIPgetOrigObjoffset(scip);
		else ret = SCIPgetTransObjoffset(scip);
	}
	return ret;
}


ObjSense PDLPModel::objSense() const
{
	SCIP_OBJSENSE sense = SCIPgetObjsense(scip);
	if (sense == -1) return ObjSense::MAX;
	else return ObjSense::MIN;
}


void PDLPModel::lbs(double* lb, int first, int last) const
{
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	std::vector<double> tempUB(ncols(), 0);
	DOMINIQS_ASSERT(SCIPlpiGetCols(lpi, first, last, lb, &tempUB[0], NULL, NULL, NULL, NULL));
}


void PDLPModel::ubs(double* ub, int first, int last) const
{
	// Done
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	std::vector<double> tempLB(ncols(), 0);
	DOMINIQS_ASSERT(SCIPlpiGetCols(lpi, 0, ncols()-1, &tempLB[0], ub, NULL, NULL, NULL, NULL));
}


void PDLPModel::objcoefs(double* obj, int first, int last) const
{
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		SCIPlpiGetObj(lpi, first, last, obj);
	}
	else
	{
		if (postsolved)
		{
			SCIP_VAR ** vars;
			int nvars = SCIPgetNOrigVars(scip);
			SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
			for ( int i = 0; i < ncols(); i++)
			{
				obj[i] = SCIPvarGetObj(vars[i]);
			}
		}
		else
		{
			SCIP_VAR ** vars;
			int nvars = SCIPgetNVars(scip);
			SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
			for ( int i = 0; i < ncols(); i++)
			{
				obj[i] = SCIPvarGetObj(vars[i]);
			}
		}
	}
}


void PDLPModel::ctypes(char* ctype, int first, int last) const
{
	// col types are obtained from the cloned problem
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	SCIP_VAR ** vars;
	int nvars;
	SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
	for (int idx = first; idx <= last; idx++)
	{
		if (SCIPvarIsIntegral(vars[idx]))
		{
			if (SCIPvarIsBinary(vars[idx])) ctype[idx] = 'B';
			else ctype[idx] = 'I';
		}
		else ctype[idx] = 'C';
	}
}


void PDLPModel::sense(char* sense, int first, int last) const
{
	// not needed
}


void PDLPModel::rhs(double* rhs, int first, int last) const
{
	// not needed
}


void PDLPModel::row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const
{
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		DOMINIQS_ASSERT((ridx >= 0) && (ridx < nrows()));
		double inft = SCIPlpiInfinity(lpi);
		double lhsval;
		double rhsval;
		int nnonz;
		std::vector<int> beg(nrows(), 0);
		std::vector<int> ind(ncols(), 0);
		std::vector<double> val(ncols(), 0);
		SCIPlpiGetRows(lpi, ridx, ridx, &lhsval, &rhsval, &nnonz , &beg[0], &ind[0], &val[0]);

		if (lhsval == -inft)
		{
			rhs = rhsval;
			sense = 'L';
		}
		else if (rhsval == inft)
		{
			rhs = lhsval;
			sense = 'G';
		}
		else if (rhsval == lhsval)
		{
			rhs = rhsval;
			sense = 'E';
		}
		else if ((rhsval < inft) && (lhsval > -inft))
		{
			rhs = rhsval;
			rngval = rhsval - lhsval;
			sense = 'R';
		}
		else
		{
			std::cout << "Something wrong happened with row: " << ridx << std::endl;
		}

		for (int idx = 0 ; idx < nnonz ; idx++)
		{
			row.push(ind[idx], val[idx]);
		}

		//  check with scip!!! (delete this when everything works)
		SCIP_CONS ** cons;
		if (isClone) cons = SCIPgetOrigConss(scip);
		else cons = SCIPgetOrigConss(scip);
		DOMINIQS_ASSERT(cons);
		unsigned int success;
		double rhsvalMIP, lhsvalMIP;
		int nnonzMIP;
		rhsvalMIP = SCIPconsGetRhs(scip, cons[ridx], &success);
		DOMINIQS_ASSERT(success);
		lhsvalMIP = SCIPconsGetLhs(scip, cons[ridx], &success);
		DOMINIQS_ASSERT(success);

		if ( sense == 'L' )
		{
			DOMINIQS_ASSERT(lhsvalMIP == -inft);
			DOMINIQS_ASSERT(rhs == rhsvalMIP);
		}
		else if ( sense == 'G' )
		{
			DOMINIQS_ASSERT(rhsvalMIP == inft);
			DOMINIQS_ASSERT(rhs == lhsvalMIP);
		}
		else if ( sense == 'E' )
		{
			DOMINIQS_ASSERT(lhsvalMIP == rhsvalMIP);
			DOMINIQS_ASSERT(rhs == lhsvalMIP);
		}
		else
		{
			DOMINIQS_ASSERT(rhs == rhsvalMIP);
			DOMINIQS_ASSERT(rngval == (rhsvalMIP - lhsvalMIP));
		}
		SCIPgetConsNVars(scip, cons[ridx], &nnonzMIP, &success);
		DOMINIQS_ASSERT((nnonz == nnonzMIP));
		DOMINIQS_ASSERT(success);
		SCIP_VAR ** vars;
		SCIPallocBufferArray(scip, &vars, nnonzMIP);
		SCIPgetConsVars(scip, cons[ridx], vars, nnonzMIP, &success);
		DOMINIQS_ASSERT(success);
		DOMINIQS_ASSERT(vars);
	}
	else
	{
		double inft = SCIPinfinity(scip);
		SCIP_CONS ** cons;
		if (isClone) cons = SCIPgetOrigConss(scip);
		else cons = SCIPgetOrigConss(scip);
		DOMINIQS_ASSERT(cons);
		SCIP_CONS * con = cons[ridx];
		unsigned int success;
		double lhsval;
		double rhsval;
		rhsval = SCIPconsGetRhs(scip, cons[ridx], &success);
		lhsval = SCIPconsGetLhs(scip, cons[ridx], &success);
		if (lhsval == -inft)
		{
			rhs = rhsval;
			sense = 'L';
		}
		else if (rhsval == inft)
		{
			rhs = lhsval;
			sense = 'G';
		}
		else if (rhsval == lhsval)
		{
			rhs = rhsval;
			sense = 'E';
		}
		else
		{
			rhs = rhsval;
			rngval = rhsval - lhsval;
			sense = 'R';
		}
		int nnzs = 0;

		SCIPgetConsNVars(scip, cons[ridx], &nnzs, &success);
		DOMINIQS_ASSERT(success);
		SCIP_VAR ** vars;
		SCIPallocBufferArray(scip, &vars, nnzs);
		SCIPgetConsVars(scip, cons[ridx], vars, nnzs, &success);
		DOMINIQS_ASSERT(success);
		DOMINIQS_ASSERT(vars);
		std::vector<int> indices(nnzs, 0);
		for (int idx = 0; idx < nnzs; idx++)
		{
			indices[idx] = SCIPvarGetProbindex(vars[idx]);
		}
		SCIPfreeBufferArray(scip, &vars);

		std::vector<double> val(ncols(), 0);
		SCIPgetConsVals(scip, cons[ridx], &val[0], nnzs, &success);
		for (int idx = 0 ; idx < nnzs ; idx++)
		{
			row.push(indices[idx], val[idx]);
		}
	}

}

void PDLPModel::rows(dominiqs::SparseMatrix& matrix) const
{
	// not needed
}


void PDLPModel::col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const
{
	// not needed
}


void PDLPModel::cols(dominiqs::SparseMatrix& matrix) const
{
	// not needed
}


void PDLPModel::colNames(std::vector<std::string>& names, int first, int last) const
{
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	names.clear();

	int nvars;
	SCIP_VAR **vars;
	if ((postsolved) || (stageLP)) SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
	else SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
	for (int idx = 0; idx < nvars; idx++)
	{
		names.push_back(SCIPvarGetName(vars[idx]));
	}


}

void PDLPModel::rowNames(std::vector<std::string>& names, int first, int last) const
{
	SCIP_CONS ** cons;
	int ncons;
	if (postsolved)
	{
		cons = SCIPgetOrigConss(scip);
		ncons = SCIPgetNOrigConss(scip);
	}
	else
	{
		cons = SCIPgetConss(scip);
		ncons = SCIPgetNConss(scip);
	}
	const char * name;
	for( int i = 0; i < ncons; i ++ )
	{
		name = SCIPconsGetName(cons[i]);
		names[i] = std::string(name);
	}
}


/* Data modifications */
void PDLPModel::addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj)
{
	// ToDo for PDLP
	double pdlpinft = solver->infinity();
	if ( ub >= 1e20 ) ub = pdlpinft;
	variables.push_back(solver->MakeNumVar( lb, ub, name ));
}


void PDLPModel::addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj)
{
	// not needed...
}


void PDLPModel::addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval)
{
	// ToDo for PDLP
	if (stageLP)
	{
		double pdlpinft = solver->infinity();
		double lhsval;
		double rhsval;
		if (sense == 'G')
		{
			rhsval = pdlpinft;
			lhsval = rhs;
		}
		else if (sense=='L')
		{
			rhsval = rhs;
			lhsval = -pdlpinft;
		}
		else if (sense=='E')
		{
			rhsval = rhs;
			lhsval = rhs;
		}
		else if (sense=='R')
		{
			rhsval = rhs;
			lhsval = rhs - rngval;
		}
		MPConstraint* const con = solver->MakeRowConstraint( lhsval, rhsval );
		for( int j = 0; j < cnt; j++ )
		{
			con->SetCoefficient( variables[idx[j]], val[j] );
		}
	}
}


void PDLPModel::delRow(int ridx)
{
	// not needed
}


void PDLPModel::delCol(int cidx)
{
	// not needed
}


void PDLPModel::delRows(int first, int last)
{
	// ToDo for PDLP
	// Just delete the whole problem and copy the lpi again...
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);
	DOMINIQS_ASSERT(solver);

	// clear solver
	solver->Clear();

	// Create model for PDLP
	double inft = SCIPlpiInfinity(lpi);
	double pdlpinft = solver->infinity();

	int nrows, ncols;

	SCIPlpiGetNCols(lpi, &ncols);
	SCIPlpiGetNRows(lpi, &nrows);

	std::vector<double> lower(ncols, 0);
	std::vector<double> upper(ncols, 0);
	std::vector<double> coefs(ncols, 0);

	objcoefs(coefs.data(), 0, -1);
	ubs(upper.data(), 0, -1);
	lbs(lower.data(), 0, -1);

	std::vector<std::string> xNames;
	colNames(xNames);

	variables = std::vector<MPVariable*>{};
	variables.reserve(ncols);

	MPObjective* const objective = solver->MutableObjective();

	for( int i = 0; i < ncols; i++ )
	{
		// ToDo check infinity values
		if (lower[i] <= -1e20) lower[i] = -pdlpinft;
		else if (upper[i] >= 1e20) upper[i] = pdlpinft;
		variables.push_back(solver->MakeNumVar( lower[i], upper[i], xNames[i] ));
	    objective->SetCoefficient( variables[i], coefs[i] );
	}

	assert( solver->NumVariables() == ncols );
	double lhsval;
	double rhsval;
	int nnonz;
	std::vector<int> beg(1, 0);
	std::vector<int> ind(ncols, 0);
	std::vector<double> val(ncols, 0);

	for( int i = 0; i < nrows; i++ )
	{
		SCIPlpiGetRows(lpi, i, i, &lhsval, &rhsval, &nnonz , &beg[0], &ind[0], &val[0]);

		MPConstraint* const con = solver->MakeRowConstraint();
		// 1e20 is infinity of SCIP
		if (lhsval <= -1e20) con->SetUB(rhsval);
		else if (rhsval >= 1e20) con->SetLB(lhsval);
		else
		{
			con->SetUB(rhsval);
			con->SetLB(lhsval);
		}
		for( int j = 0; j < nnonz; j++ )
		{
			con->SetCoefficient(variables[ind[j]], val[j]);
		}
	}
    assert( solver->NumConstraints() == nrows);
    objective->SetMinimization();
}


void PDLPModel::delCols(int first, int last)
{
	// ToDo for PDLP
}


void PDLPModel::objSense(ObjSense objsen)
{
	// ToDo remove lpi
	MPObjective* const objective = solver->MutableObjective();
	if (objsen == ObjSense::MIN)
	{
		SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE);
		SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE);
    	objective->SetMinimization();
	}
	else
	{
		SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MAXIMIZE);
		SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);
    	objective->SetMaximization();
	}

}


void PDLPModel::objOffset(double val)
{
	// Done
	// Add offset to MIP
	SCIPaddOrigObjoffset(scip, val);
}

void PDLPModel::lb(int cidx, double val)
{
	SCIP_VAR **vars;
	SCIPgetVarsData(scip, &vars, NULL, NULL, NULL, NULL, NULL);
	SCIPchgVarLb(scip,vars[cidx],val);
}


void PDLPModel::lbs(int cnt, const int* cols, const double* values)
{
	// not needed
}


void PDLPModel::ub(int cidx, double val)
{
	SCIP_VAR **vars;
	SCIPgetVarsData(scip, &vars, NULL, NULL, NULL, NULL, NULL);
	SCIPchgVarUb(scip,vars[cidx],val);
	// ToDo for PDLP (needed onl when we use stage 3)
}


void PDLPModel::ubs(int cnt, const int* cols, const double* values)
{
	// not needed
}


void PDLPModel::fixCol(int cidx, double val)
{
	// not needed
}


void PDLPModel::objcoef(int cidx, double val)
{
	// not needed
}


void PDLPModel::objcoefs(int cnt, const int* cols, const double* values)
{
	if (stageLP)
	{
		MPObjective* const objective = solver->MutableObjective();
		for( int i = 0; i < ncols(); i++ )
		{
			objective->SetCoefficient( variables[i], values[i] );
		}
	}

}


void PDLPModel::ctype(int cidx, char val)
{
	// not needed...
	// change column type (not needed since we created the MIP and LP separately.)
}


void PDLPModel::ctypes(int cnt, const int* cols, const char* values)
{
	// not needed...
	// change column type for multiple columns (not needed since we created the MIP and LP separately.)
}


void PDLPModel::switchToLP()
{
	// Done
	stageLP = true;
}

void PDLPModel::switchToMIP()
{
	// Done
	stageLP = false;
	postsolved = false;
	nvarsStartStage3 = ncols();
	nconsStartStage3 = nrows();
}

/* Private interface */

PDLPModel* PDLPModel::clone_impl() const
{
	DOMINIQS_ASSERT(scip);
	// solve the root node to get the LP (only root, only one lp iteration, no presolve)
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "limits/nodes", 1));
	DOMINIQS_ASSERT(SCIPsetRealParam(scip, "limits/gap", SCIPinfinity(scip)));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "lp/iterlim", 1));
	DOMINIQS_ASSERT(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
	DOMINIQS_ASSERT((SCIPsolve(scip) == SCIP_OKAY));
	DOMINIQS_ASSERT(SCIPresetParams(scip));
	double offset = objOffset();
	// create copy
	PDLPModel* copy = new PDLPModel();
	DOMINIQS_ASSERT(SCIPcreate(&copy->scip));
	DOMINIQS_ASSERT(SCIPcopy(scip, copy->scip, NULL, NULL, NULL, true, false, false, false, NULL));
	DOMINIQS_ASSERT(copy->scip);
	copy->isClone = true;
	DOMINIQS_ASSERT(copy->lpi);
	// ToDo implementation for PDLP is missing. clone_impl() is only called when presolving is disabled.
	return copy;

}
PDLPModel* PDLPModel::presolvedmodel_impl()
{
	DOMINIQS_ASSERT(scip);
	// solve the root node to get the LP (only root, only one lp iteration, no presolve since it is already presolved)
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
	DOMINIQS_ASSERT(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "limits/nodes", 1));
	DOMINIQS_ASSERT(SCIPsetRealParam(scip, "limits/gap", SCIPinfinity(scip)));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "lp/iterlim", 1));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "lp/rootiterlim", 1));

	dominiqs::StopWatch chrono;
	dominiqs::StopWatch chronoRows;
	dominiqs::StopWatch chronoCols;

	chrono.start();
	
	std::cout << "extracting the LP" << std::endl;
	SCIPsolve(scip);
	DOMINIQS_ASSERT(SCIPresetParams(scip));
	double offset = objOffset();
	// create copy
	PDLPModel* copy = new PDLPModel();
	DOMINIQS_ASSERT(SCIPcreate(&copy->scip));
	DOMINIQS_ASSERT(SCIPcopy(scip, copy->scip, NULL, NULL, NULL, true, false, false, false, NULL));
	DOMINIQS_ASSERT(copy->scip);
	DOMINIQS_ASSERT(SCIPaddOrigObjoffset(copy->scip, offset));
	// get LPI
	DOMINIQS_ASSERT(SCIPgetLPI(scip, &lpi));
	SCIP_OBJSEN senseLP;
	DOMINIQS_ASSERT(SCIPlpiGetObjsen(lpi, &senseLP));
	// create copylpi
	DOMINIQS_ASSERT(SCIPlpiCreate(&copy->lpi, NULL, "copylpi", senseLP));
	copyLpi(&copy->lpi, lpi);
	copy->isClone = true;
	DOMINIQS_ASSERT(copy->lpi);
	std::cout << "LP is extracted" << std::endl;

	std::cout << "Creating PDLP model" << std::endl;

	double pdlpinft = std::numeric_limits<double>::infinity();
	double scipinf = 1e20;

	// Create model for PDLP
	// ToDo write this as a function
	SCIPlpiGetNCols(lpi, &numColsPresolved);
	SCIPlpiGetNRows(lpi, &numRowsPresolved);

	std::vector<double> lower(numColsPresolved, 0);
	std::vector<double> upper(numColsPresolved, 0);
	std::vector<double> coefs(numColsPresolved, 0);

	objcoefs(coefs.data(), 0, -1);
	ubs(upper.data(), 0, -1);
	lbs(lower.data(), 0, -1);

	std::vector<std::string> xNames;
	colNames(xNames);

	variables = std::vector<MPVariable*>{};
	variables.reserve(numColsPresolved);

	MPObjective* const objective = solver->MutableObjective();

	chronoCols.start();
	for( int i = 0; i < numColsPresolved; i++ )
	{
		if (lower[i] <= -scipinf) lower[i] = -pdlpinft;
		if (upper[i] >= scipinf) upper[i] = pdlpinft;
		variables.push_back(solver->MakeNumVar( lower[i], upper[i], xNames[i] ));
	    objective->SetCoefficient( variables[i], coefs[i] );
	}
	chronoCols.stop();
	std::cout << "PDLP adding columns in model took " << chronoCols.getElapsed() << " seconds" << std::endl;

	assert( solver->NumVariables() == numColsPresolved );

	chronoRows.start();

	for( int i = 0; i < numRowsPresolved; i++ )
	{
		double lhsval;
		double rhsval;
		int nnonz;
		std::vector<int> beg(1, 0);
		std::vector<int> ind(numColsPresolved, 0);
		std::vector<double> val(numColsPresolved, 0);

		SCIPlpiGetRows(lpi, i, i, &lhsval, &rhsval, &nnonz , &beg[0], &ind[0], &val[0]);
		MPConstraint* const con = solver->MakeRowConstraint();
		if (lhsval <= -scipinf) con->SetUB(rhsval);
		else if (rhsval >= scipinf) con->SetLB(lhsval);
		else
		{
			con->SetUB(rhsval);
			con->SetLB(lhsval);
		}
		for( int j = 0; j < nnonz; j++ )
		{
			con->SetCoefficient(variables[ind[j]], val[j]);
		}
	}
    assert(solver->NumConstraints() == numRowsPresolved);
	chronoRows.stop();
	std::cout << "PDLP adding rows in model took " << chronoRows.getElapsed() << " seconds" << std::endl;
	
    objective->SetMinimization();
	copy->numColsPresolved = numColsPresolved;
	copy->numRowsPresolved = numRowsPresolved;
	copy->verbosity = verbosity;

	copy->solver = solver;
	copy->variables = variables;

	// ToDo should switch to qp and not the proto model to avoid the overhead of copying the model in each iteration

	params = operations_research::pdlp::PrimalDualHybridGradientParams();

	copy->params = params;
	// Default values for pdlp
	copy->pdlpTol = 1.0e-6;
	copy->pdlpTolDecreaseFactor = 1.0;
 	copy->pdlpWarmStart = 0;

	chrono.stop();
	std::cout << "PDLP model created in " << chrono.getElapsed() << " seconds" << std::endl;

	return copy;
}


void PDLPModel::copyLpi(SCIP_LPI** targetlpi, SCIP_LPI* srclpi)
{
    // Done
	std::vector<double> val, lb, ub , lhs, rhs, obj;
	std::vector<int> beg, ind;
	int nnonz, nrows, ncols;
	SCIP_OBJSEN senseLP;

	SCIPlpiGetObjsen(srclpi, &senseLP);
	SCIPlpiGetNNonz(srclpi, &nnonz);
	SCIPlpiGetNRows(srclpi, &nrows);
	SCIPlpiGetNCols(srclpi, &ncols);

	val.resize(nnonz);
	ind.resize(nnonz);
	beg.resize(ncols);
	lb.resize(ncols);
	ub.resize(ncols);
	lhs.resize(nrows);
	rhs.resize(nrows);
	obj.resize(ncols);

   	SCIPlpiGetCols(srclpi, 0, ncols - 1, lb.data(), ub.data(), &nnonz, beg.data(), ind.data(), val.data());
	SCIPlpiGetSides(srclpi, 0, nrows-1, lhs.data(), rhs.data());
	SCIPlpiGetObj(srclpi, 0, ncols-1, obj.data());
	DOMINIQS_ASSERT(SCIPlpiLoadColLP(*targetlpi, senseLP, ncols, obj.data(), lb.data(), ub.data(), NULL, nrows, lhs.data(), rhs.data(), NULL, nnonz, beg.data(), ind.data(), val.data()));

}
