/**
 * @file mipmodel.h
 * @brief MIP Model Interface
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2019
 */

#ifndef MIPMODEL_H
#define MIPMODEL_H

#include <string>
#include <vector>
#include <memory>
#include <utils/maths.h>


/* Objective Sense values */
enum class ObjSense {
	MIN = 1,
	MAX = -1
};


enum class IntParam {
	Threads,
	SolutionLimit,
	NodeLimit,
	IterLimit
};


enum class DblParam {
	TimeLimit,
	FeasibilityTolerance,
	IntegralityTolerance
};


enum class IntAttr {
	Nodes,
	NodesLeft,
	BarrierIterations,
	SimplexIterations
};


enum class DblAttr {
	MIPDualBound
};


/* Interface for a MIP model and solver */
class MIPModelI
{
public:
	virtual ~MIPModelI() {}
	std::unique_ptr<MIPModelI> clone() const { return std::unique_ptr<MIPModelI>(this->clone_impl()); }
	/* Read/Write */
	virtual void readModel(const std::string& filename) = 0;
	virtual void writeModel(const std::string& filename, const std::string& format="") const = 0;
	virtual void writeSol(const std::string& filename) const = 0;
	/* Solve */
	virtual void lpopt(char method) = 0;
	virtual void mipopt() = 0;
	/* Presolve/Postsolve */
	virtual void presolve() = 0;
	virtual void postsolve() = 0;
	std::unique_ptr<MIPModelI> presolvedModel() { return std::unique_ptr<MIPModelI>(this->presolvedmodel_impl()); }
	virtual std::vector<double> postsolveSolution(const std::vector<double>& preX) const = 0;
	/* Get solution */
	virtual double objval() const = 0;
	virtual void sol(double* x, int first = 0, int last = -1) const = 0;
	virtual bool isPrimalFeas() const = 0;
	/* Parameters */
	virtual void handleCtrlC(bool flag) = 0;
	virtual bool aborted() const = 0;
	virtual void seed(int seed) = 0;
	virtual void logging(bool log) = 0;
	virtual int intParam(IntParam which) const = 0;
	virtual void intParam(IntParam which, int value) = 0;
	virtual double dblParam(DblParam which) const = 0;
	virtual void dblParam(DblParam which, double value) = 0;
	virtual int intAttr(IntAttr which) const = 0;
	virtual double dblAttr(DblAttr which) const = 0;
	/* Access model data */
	virtual int nrows() const = 0;
	virtual int ncols() const = 0;
	virtual int nnz() const = 0;
	virtual double objOffset() const = 0;
	virtual ObjSense objSense() const = 0;
	virtual void lbs(double* lb, int first = 0, int last = -1) const = 0;
	virtual void ubs(double* ub, int first = 0, int last = -1) const = 0;
	virtual void objcoefs(double* obj, int first = 0, int last = -1) const = 0;
	virtual void ctypes(char* ctype, int first = 0, int last = -1) const = 0;
	virtual void sense(char* sense, int first = 0, int last = -1) const = 0;
	virtual void rhs(double* rhs, int first = 0, int last = -1) const = 0;
	virtual void row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const = 0;
	virtual void rows(dominiqs::SparseMatrix& matrix) const = 0;
	virtual void col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const = 0;
	virtual void cols(dominiqs::SparseMatrix& matrix) const = 0;
	virtual void colNames(std::vector<std::string>& names, int first = 0, int last = -1) const = 0;
	virtual void rowNames(std::vector<std::string>& names, int first = 0, int last = -1) const = 0;
	/* Data modifications */
	virtual void addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj) = 0;
	virtual void addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj) = 0;
	virtual void addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval = 0.0) = 0;
	virtual void delRow(int ridx) = 0;
	virtual void delCol(int cidx) = 0;
	virtual void delRows(int first, int last) = 0;
	virtual void delCols(int first, int last) = 0;
	virtual void objSense(ObjSense objsen) = 0;
	virtual void objOffset(double val) = 0;
	virtual void lb(int cidx, double val) = 0;
	virtual void lbs(int cnt, const int* cols, const double* values) = 0;
	virtual void ub(int cidx, double val) = 0;
	virtual void ubs(int cnt, const int* cols, const double* values) = 0;
	virtual void fixCol(int cidx, double val) = 0;
	virtual void objcoef(int cidx, double val) = 0;
	virtual void objcoefs(int cnt, const int* cols, const double* values) = 0;
	virtual void ctype(int cidx, char val) = 0;
	virtual void ctypes(int cnt, const int* cols, const char* values) = 0;
	virtual void switchToLP() = 0;
private:
	virtual MIPModelI* clone_impl() const = 0;
	virtual MIPModelI* presolvedmodel_impl() = 0;
};


using MIPModelPtr = std::shared_ptr<MIPModelI>;


#endif /* MIPMODEL_H */
