#include "feaspump/scipmodel.h"
#include <signal.h>
#include <cstring>
#include <climits>
#include <numeric>
#include <fmt/format.h>
#include <cassert>
#include <utils/maths.h>


SCIP_STAGE SCIPModel::getProbStage()
{
	DOMINIQS_ASSERT(scip);
	return SCIPgetStage(scip);
}


std::vector<std::string> SCIPModel::getVarNames()
{
	DOMINIQS_ASSERT(scip);
	SCIP_STAGE stage = SCIPgetStage(scip);
	if ( stage == SCIP_STAGE_PROBLEM) ;
	{
		colNames(origColNames);
		return origColNames;
	}	
	if ( stage == SCIP_STAGE_TRANSFORMED) ;
	{
		colNames(transColNames);
		return transColNames;
	}		
	if ( stage == SCIP_STAGE_PRESOLVED) ;
	{
		colNames(presColNames);
		return presColNames;
	}			
}

SCIPModel::SCIPModel()
{	
	// Done
	// initialize scip problem
	DOMINIQS_ASSERT(SCIPcreate(&scip));
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_PROBLEM));
	DOMINIQS_ASSERT(SCIPincludeDefaultPlugins(scip));
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
}


SCIPModel::SCIPModel(SCIP* _prob, bool _ownProb) : scip(_prob), ownProb(_ownProb)
{	
	// Done
	// assert problem
	DOMINIQS_ASSERT(scip);
}

/* Destroyer */
SCIPModel::~SCIPModel()
{	
	// // ToDo
	// DOMINIQS_ASSERT(scip);
	// SCIPfreeTransform(scip);
	// SCIPfree(&scip);
	// if (stageLP)
	// {
	// 	DOMINIQS_ASSERT(lpi);
	// 	SCIPlpiFree(&lpi);
	// }
}


/* Read/Write */
void SCIPModel::readModel(const std::string& filename)
{	
	// Done
	// read problem from file
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(SCIPreadProb(scip, filename.c_str(), NULL));
	// colNames(origColNames);
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_PROBLEM));
	DOMINIQS_ASSERT(SCIPtransformProb(scip));
	DOMINIQS_ASSERT((getProbStage(),SCIP_STAGE_TRANSFORMED));
	// colNames(transColNames);

}

void SCIPModel::writeModel(const std::string& filename, const std::string& format) const
{
    // later
}

void SCIPModel::writeSol(const std::string& filename) const
{
    // later
}


/* Solve */
double SCIPModel::lpopt(char method, bool decrease_tol, bool initial)
{	
	// Done
	// Which is the default method for soplex?
	DOMINIQS_ASSERT(lpi);
	switch(method)
	{
		// should the dual simplex be the default method?
		case 'S': SCIPlpiSolveDual(lpi); break;

		case 'P': SCIPlpiSolvePrimal(lpi); break;
		case 'D': SCIPlpiSolveDual(lpi); break;
		case 'B': SCIPlpiSolveBarrier(lpi, 1); break;
		default: throw std::runtime_error("Unexpected method for lpopt");
	}
	return 0.0;

	// // some checks
	// int solved;
	// int opt;
	// opt = SCIPlpiIsOptimal(lpi);
	// solved = SCIPlpiWasSolved(lpi);
	// DOMINIQS_ASSERT(opt);
	// DOMINIQS_ASSERT(solved);
}

/* Solve mip */
void SCIPModel::mipopt()
{	
	// Done
	// solve mip
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(SCIPsolve(scip));
}

/* presolve mip */
void SCIPModel::presolve()
{	
	// presolve the problem
	DOMINIQS_ASSERT(scip);
	// SCIPsetIntParam(scip, "display/verblevel", 4);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/indicator" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/knapsack" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/logicor" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/setppc" , false);	
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/varbound" , false);
	SCIPsetBoolParam(scip,"constraints/linear/upgrade/xor" , false);
	DOMINIQS_ASSERT(SCIPpresolve(scip));
	DOMINIQS_ASSERT(SCIPresetParams(scip));
	SCIPsetIntParam(scip, "display/verblevel", 0);
}

// Postsolve the problem
void SCIPModel::postsolve()
{	
	// Done
	// probably not needed since we can access original vars and constraints in scip
	if (!stageLP) 
	{
		postsolved = true;
		DOMINIQS_ASSERT(scip);
	}
}

/* get solution vector in the original space */
std::vector<double> SCIPModel::postsolveSolution(const std::vector<double>& preX) const
{	
	// Done
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
double SCIPModel::objval() const
{
	DOMINIQS_ASSERT(lpi);
	double objval;
	SCIPlpiGetObjval(lpi, &objval);
	return objval;
}


void SCIPModel::sol(double* x, int first, int last) const
{	
	// Done (check stage 3 again)

	if (stageLP)
	{
		DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
		if (last == -1)  last = ncols()-1;
		DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
		DOMINIQS_ASSERT(lpi);
		std::vector<double> projSol;
		projSol.resize(ncols(),0);
		SCIPlpiGetSol(lpi, NULL, &projSol[0], NULL, NULL, NULL);
		for (int idx = 0; idx <= last; idx++) x[idx] = projSol[idx];
	}
	else
	{
		DOMINIQS_ASSERT(scip);
		// get solution
		SCIP_SOL *sol;
		sol = SCIPgetBestSol(scip);
		// retransform to original solution space
		SCIPretransformSol(scip, sol);
		// get original variables
		SCIP_VAR **vars;
		int nvars = SCIPgetNOrigVars(scip);
		SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
		// get solution vector x
		for ( int i = first; i <= last ; i++)
		{	
			x[i] = SCIPgetSolVal(scip, sol, vars[i]);
		}
		SCIPfreeSol(scip, &sol);
	}


}


bool SCIPModel::isPrimalFeas() const
{	
	// Done
	if (stageLP)
	{	
		DOMINIQS_ASSERT(lpi);
		return (SCIPlpiWasSolved(lpi) && (SCIPlpiIsOptimal(lpi) > 0));
	}
	else
	{
		SCIP_SOL * sol = SCIPgetBestSol(scip);
		if ( sol == NULL) return 0;
		else return 1;
	}

}



/* Parameters */
void SCIPModel::handleCtrlC(bool flag)
{
	// ToDo
	// SCIPpressedCtrlC
}


bool SCIPModel::aborted() const
{
	// ToDo
    return false;
}


void SCIPModel::seed(int seed)
{
	// Done
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);
	SCIPsetIntParam(scip,"branching/random/seed", seed);
	SCIPlpiSetIntpar(lpi, SCIP_LPPAR_RANDOMSEED, seed);
}


void SCIPModel::logging(bool log)
{
	// Done (only for scip)
	DOMINIQS_ASSERT(scip);
	if (log)  SCIPsetIntParam(scip, "display/verblevel", 4);
	else      SCIPsetIntParam(scip, "display/verblevel", 0);
}

int SCIPModel::intParam(IntParam which) const
{
	// ToDo
	
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


void SCIPModel::intParam(IntParam which, int value)
{
	// ToDo
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
				SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, value);
				break;
			case IntParam::PdlpWarmStart:
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
			default:
				throw std::runtime_error("Unknown integer parameter");
		}
	}

}


double SCIPModel::dblParam(DblParam which) const
{	
	// ToDo
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
			default:
				throw std::runtime_error("Unknown double parameter");
		}
	}

	return value;
}


void SCIPModel::dblParam(DblParam which, double value)
{
	// ToDo
	if (stageLP)
	{
		switch(which)
		{
			case DblParam::TimeLimit:
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LPTILIM, value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, value);
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
				SCIPsetRealParam(scip, "limits/time", value);
				break;
			case DblParam::FeasibilityTolerance:
				SCIPsetRealParam(scip, "numerics/feastol", value);
				break;
			case DblParam::IntegralityTolerance:
				SCIPsetRealParam(scip, "numerics/feastol", value);
				break;
			default:
				throw std::runtime_error("Unknown double parameter");
		}

	}
	
}


int SCIPModel::intAttr(IntAttr which) const
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
				break;
			case IntAttr::PDLPIterations:
				break;
			default:
				throw std::runtime_error("Unknown integer attribute");
		}
	}
    return value;
}


double SCIPModel::dblAttr(DblAttr which) const
{
	// later
    return 0.0;
}

void SCIPModel::terminationReason(std::string& reason)
{
	// ToDo
	reason = "-";
}

/* Access model data */
int SCIPModel::nrows() const
{
	// Done
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
			DOMINIQS_ASSERT(SCIPlpiGetNRows(lpi, &ret));
		}
		else 
		{
			if (postsolved) ret = nconsStartStage3;
			else ret = SCIPgetNConss(scip);		
		}
	}
	return ret;
}


int SCIPModel::ncols() const
{
	// Done
    int ret;
	if (!isClone) 
	{	
		if (postsolved) ret = SCIPgetNOrigVars(scip);
		else ret = SCIPgetNVars(scip);
	}
	else
	{
		if (stageLP)
		{
			DOMINIQS_ASSERT(SCIPlpiGetNCols(lpi, &ret));
		} 
		else 
		{
			if (postsolved) ret = nvarsStartStage3;
			else ret = SCIPgetNVars(scip);		
		}
	}
	return ret;
}


int SCIPModel::nnz() const
{
	// ToDo
    int ret;
	if (stageLP) DOMINIQS_ASSERT(SCIPlpiGetNNonz(lpi, &ret));
	else ret = SCIPgetNNZs(scip);
	return ret;
}

double SCIPModel::objOffset() const
{
	// Done
	double ret = 0.0;
	if (isClone) ret = SCIPgetOrigObjoffset(scip);
	else 
	{
		if (postsolved) ret = SCIPgetOrigObjoffset(scip);
		else ret = SCIPgetTransObjoffset(scip);
	}
	return ret;
}


ObjSense SCIPModel::objSense() const
{
	// Done
	SCIP_OBJSENSE sense = SCIPgetObjsense(scip);
	if (sense == -1) return ObjSense::MAX;
	else return ObjSense::MIN;
}


void SCIPModel::lbs(double* lb, int first, int last) const
{	
	// Done
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	std::vector<double> tempUB(ncols(), 0);
	DOMINIQS_ASSERT(SCIPlpiGetCols(lpi, first, last, lb, &tempUB[0], NULL, NULL, NULL, NULL));

	// // double-check with scip 
	// SCIP_VAR ** vars;
	// int nvars = SCIPgetNOrigVars(scip);
	// int nbinvars;
	// int nintvars;
	// int nimplvars;
	// int ncontvars;
	// SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, &ncontvars);
	// for (int i=0; i < nvars; i++)
	// {
	// 	DOMINIQS_ASSERT((lb[i] == SCIPvarGetLbOriginal(vars[i])));
	// }
}


void SCIPModel::ubs(double* ub, int first, int last) const
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
	
	// // double-check with scip 
	// SCIP_VAR ** vars;
	// int nvars = SCIPgetNOrigVars(scip);
	// int nbinvars;
	// int nintvars;
	// int nimplvars;
	// int ncontvars;
	// SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, &ncontvars);
	// for (int i=0; i < nvars; i++)
	// {
	// 	std::cout << i << " " << ub[i] << " " << SCIPvarGetUbOriginal(vars[i]) << std::endl;
	// 	DOMINIQS_ASSERT((ub[i] == SCIPvarGetUbOriginal(vars[i])));
	// }
}


void SCIPModel::objcoefs(double* obj, int first, int last) const
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
		// check with scip!!!
		// SCIP_VAR ** vars;
		// int nvars;
		// SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
		// for ( int i = 0; i < ncols(); i++)
		// {	
		// 	DOMINIQS_ASSERT(( obj[i] == SCIPvarGetObj(vars[i]) ));
		// }

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


void SCIPModel::ctypes(char* ctype, int first, int last) const
{
	// Done
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


void SCIPModel::sense(char* sense, int first, int last) const
{
	// later
}


void SCIPModel::rhs(double* rhs, int first, int last) const
{
	// later
}


void SCIPModel::row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const
{
	// ToDo
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


void SCIPModel::rows(dominiqs::SparseMatrix& matrix) const
{
	// later
}


void SCIPModel::col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const
{
	// later
}


void SCIPModel::cols(dominiqs::SparseMatrix& matrix) const
{
	// later
}


void SCIPModel::colNames(std::vector<std::string>& names, int first, int last) const
{
	// Done
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

void SCIPModel::rowNames(std::vector<std::string>& names, int first, int last) const
{
	// Done
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
void SCIPModel::addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj)
{
	// Done
	DOMINIQS_ASSERT(scip);
	DOMINIQS_ASSERT(lpi);

	char* cname = (char*)(name.c_str());

	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		// ToDo add empty column to lp
		SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, &cname, 0, NULL, NULL, NULL);

	}
	else
	{
		SCIP_VAR* var;
		SCIPcreateVarBasic(scip,&var, cname, lb, ub, obj, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip, var);
		SCIPreleaseVar(scip, &var);
	}

}


void SCIPModel::addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj)
{
	// later...
}


void SCIPModel::addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval)
{

	// Done
	if (stageLP)
	{
		double inft = SCIPlpiInfinity(lpi);
		double lhsval;
		double rhsval;
		if (sense == 'G')
		{
			rhsval = inft;
			lhsval = rhs;
		}
		else if (sense=='L')
		{
			rhsval = rhs;
			lhsval = -inft;
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
		char* rname = (char*)(name.c_str());
		int nnonz = 2;
		int beg = 0;

		int ind;
		double vall;
		int nnnrows = nrows();
		SCIPlpiAddRows(lpi, 1, &lhsval, &rhsval, &rname, nnonz, &beg, idx, val);
	}
	else
	{
		double inft = SCIPinfinity(scip);
		double lhsval;
		double rhsval;
		if (sense == 'G')
		{
			rhsval = inft;
			lhsval = rhs;
		}
		else if (sense=='L')
		{
			rhsval = rhs;
			lhsval = -inft;
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

		SCIP_VAR **vars;
		int nnonz = 2;
		int nvars;
		SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);


		SCIP_Real consvals[] = {val[0], val[1]};
		SCIP_VAR** consvars;
		SCIPallocBufferArray(scip, &consvars, 2);

		consvars[0] = vars[idx[0]];
		consvars[1] = vars[idx[1]];

		char* rname = (char*)(name.c_str());
		SCIP_CONS *cons;
		SCIPcreateConsBasicLinear(scip, &cons, rname, 2 ,consvars, consvals, lhsval, rhsval);
		SCIPaddCons(scip, cons);
		SCIPreleaseCons(scip, &cons);
		SCIPfreeBufferArray(scip, &consvars);
	}
}


void SCIPModel::delRow(int ridx)
{
	// later
	DOMINIQS_ASSERT(lpi);
	SCIPlpiDelRows(lpi, ridx, ridx);
}


void SCIPModel::delCol(int cidx)
{
	// later
	DOMINIQS_ASSERT(lpi);
	SCIPlpiDelCols(lpi, cidx, cidx);
}


void SCIPModel::delRows(int first, int last)
{
	// Done (not for stage 3...)
	if (stageLP) 
	{
		DOMINIQS_ASSERT(lpi);
		DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
		DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
		DOMINIQS_ASSERT(first <= last);
		SCIPlpiDelRows(lpi, first, last);
	}
	else
	{
		// SCIP_CONS ** cons;
		// cons = SCIPgetConss(scip);
		// for ( int i = 0; i <= last; i++)
		// {
		// 	SCIPdelCons(scip, cons[i]);
		// }
	}

}


void SCIPModel::delCols(int first, int last)
{	
	// Done (not possible for stage 3...)
	if (stageLP)
	{
		DOMINIQS_ASSERT(lpi);
		DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
		DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
		DOMINIQS_ASSERT(first <= last);
		SCIPlpiDelCols(lpi, first, last);
	} 
	else
	{
		// SCIP_VAR **vars;
		// int nvars;
		// SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
		// SCIP_Bool deleted;
		// for (int idx = 0; idx <= last; idx++)
		// {
		// 	DOMINIQS_ASSERT(SCIPdelVar(scip, vars[idx], &deleted));
		// 	DOMINIQS_ASSERT(deleted);
		// }
	}
}


void SCIPModel::objSense(ObjSense objsen)
{
	// ToDo change objsense for Stage3 too!
	if (objsen == ObjSense::MIN) 
	{
		SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE);
		SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE);
	}
	else 
	{
		SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MAXIMIZE);
		SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);

	}
	
}


void SCIPModel::objOffset(double val)
{	
	// Done
	// Add offset to MIP
	SCIPaddOrigObjoffset(scip, val);
}

void SCIPModel::lb(int cidx, double val)
{
	SCIP_VAR **vars;
	SCIPgetVarsData(scip, &vars, NULL, NULL, NULL, NULL, NULL);
	SCIPchgVarLb(scip,vars[cidx],val);
}


void SCIPModel::lbs(int cnt, const int* cols, const double* values)
{
	// later
}


void SCIPModel::ub(int cidx, double val)
{
	SCIP_VAR **vars;
	SCIPgetVarsData(scip, &vars, NULL, NULL, NULL, NULL, NULL);
	SCIPchgVarUb(scip,vars[cidx],val);
}


void SCIPModel::ubs(int cnt, const int* cols, const double* values)
{
	// later
}


void SCIPModel::fixCol(int cidx, double val)
{
	// later
}


void SCIPModel::objcoef(int cidx, double val)
{
	// later
}


void SCIPModel::objcoefs(int cnt, const int* cols, const double* values)
{
	// set new objective
	if (stageLP) SCIPlpiChgObj(lpi, ncols(), cols, values);
	else
	{
		if (!postsolved)
		{
			SCIP_VAR **vars;
			int nvars = SCIPgetNVars(scip);
			SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL);
			// set objective values for the variables
			for (int idx = 0; idx < cnt; idx++)
			{
				SCIPchgVarObj(scip, vars[idx], values[idx]);
			}
		}
	}

}


void SCIPModel::ctype(int cidx, char val)
{
	// later...
	// change column type (not needed since we created the MIP and LP separately.)
}


void SCIPModel::ctypes(int cnt, const int* cols, const char* values)
{
	// later...
	// change column type for multiple columns (not needed since we created the MIP and LP separately.)

}


void SCIPModel::switchToLP()
{	
	// Done
	stageLP = true;
}

void SCIPModel::switchToMIP()
{
	// Done
	stageLP = false;
	postsolved = false;
	nvarsStartStage3 = ncols();
	nconsStartStage3 = nrows();

}

/* Private interface */

SCIPModel* SCIPModel::clone_impl() const
{
    // Done (but hacky!)
	DOMINIQS_ASSERT(scip);
	// solve the root node to get the LP (only root, only one lp iteration, no presolve)
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "limits/nodes", 1));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "lp/iterlim", 1));
	DOMINIQS_ASSERT(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
	DOMINIQS_ASSERT((SCIPsolve(scip) == SCIP_OKAY));
	DOMINIQS_ASSERT(SCIPresetParams(scip));
	double offset = objOffset();
	// create copy
	SCIPModel* copy = new SCIPModel();
	DOMINIQS_ASSERT(SCIPcreate(&copy->scip));
	DOMINIQS_ASSERT(SCIPcopy(scip, copy->scip, NULL, NULL, NULL, true, false, false, false, NULL));
	DOMINIQS_ASSERT(copy->scip);
	DOMINIQS_ASSERT(SCIPaddOrigObjoffset(copy->scip, offset));
	// get LPI for clone
	DOMINIQS_ASSERT(SCIPgetLPI(scip, &lpi));
	DOMINIQS_ASSERT(SCIPgetLPI(scip, &copy->lpi));
	copy->isClone = true;
	DOMINIQS_ASSERT(copy->lpi);
	return copy;

}
SCIPModel* SCIPModel::presolvedmodel_impl()
{
    // Done (but hacky!)
	DOMINIQS_ASSERT(scip);
	// solve the root node to get the LP (only root, only one lp iteration, no presolve since it is already presolved)
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "limits/nodes", 1));
	DOMINIQS_ASSERT(SCIPsetLongintParam(scip, "lp/iterlim", 0));
	DOMINIQS_ASSERT(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
	DOMINIQS_ASSERT(SCIPreadParams(scip, "FPscip.set"));
	DOMINIQS_ASSERT((SCIPsolve(scip) == SCIP_OKAY));
	DOMINIQS_ASSERT(SCIPresetParams(scip));
	double offset = objOffset();
	// create copy
	SCIPModel* copy = new SCIPModel();
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
	return copy;
}


void SCIPModel::copyLpi(SCIP_LPI** targetlpi, SCIP_LPI* srclpi)
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

SCIP_VAR ** SCIPModel::varsData()
{
	return varsAtStart;
}

