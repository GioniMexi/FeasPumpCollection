/**
 * @file xprsmodel.h
 * @brief Implementation of MIPModelI for XPRESS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2019
 */

#ifndef XPRSMODEL_H
#define XPRSMODEL_H

#include "mipmodel.h"
#include <xprs.h>

class XPRSModel : public MIPModelI
{
public:
	XPRSModel();
	XPRSModel(XPRSprob prob, bool _ownProb = false);
	~XPRSModel() override;
	std::unique_ptr<XPRSModel> clone() const { return std::unique_ptr<XPRSModel>(this->clone_impl()); }
	/* Read/Write */
	void readModel(const std::string& filename) override;
	void writeModel(const std::string& filename, const std::string& format="") const override;
	void writeSol(const std::string& filename) const override;
	/* Solve */
	void lpopt(char method) override;
	void mipopt() override;
	/* Presolve/postsolve */
	void presolve() override;
	void postsolve() override;
	std::unique_ptr<XPRSModel> presolvedModel() { return std::unique_ptr<XPRSModel>(this->presolvedmodel_impl()); }
	std::vector<double> postsolveSolution(const std::vector<double>& preX) const override;
	/* Get solution */
	double objval() const override;
	void sol(double* x, int first = 0, int last = -1) const override;
	bool isPrimalFeas() const override;
	/* Parameters */
	void handleCtrlC(bool flag) override;
	bool aborted() const override;
	void seed(int seed) override;
	void logging(bool log) override;
	int intParam(IntParam which) const override;
	void intParam(IntParam which, int value) override;
	double dblParam(DblParam which) const override;
	void dblParam(DblParam which, double value) override;
	int intAttr(IntAttr which) const override;
	double dblAttr(DblAttr which) const override;
	/* Access model data */
	int nrows() const override;
	int ncols() const override;
	int nnz() const override;
	double objOffset() const override;
	ObjSense objSense() const override;
	void lbs(double* lb, int first = 0, int last = -1) const override;
	void ubs(double* ub, int first = 0, int last = -1) const override;
	void objcoefs(double* obj, int first = 0, int last = -1) const override;
	void ctypes(char* ctype, int first = 0, int last = -1) const override;
	void sense(char* sense, int first = 0, int last = -1) const override;
	void rhs(double* rhs, int first = 0, int last = -1) const override;
	void row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const override;
	void rows(dominiqs::SparseMatrix& matrix) const override;
	void col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const override;
	void cols(dominiqs::SparseMatrix& matrix) const override;
	void colNames(std::vector<std::string>& names, int first = 0, int last = -1) const override;
	void rowNames(std::vector<std::string>& names, int first = 0, int last = -1) const override;
	/* Data modifications */
	void addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj) override;
	void addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj) override;
	void addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval = 0.0) override;
	void delRow(int ridx) override;
	void delCol(int cidx) override;
	void delRows(int first, int last) override;
	void delCols(int first, int last) override;
	void objSense(ObjSense objsen) override;
	void objOffset(double val) override;
	void lb(int cidx, double val) override;
	void lbs(int cnt, const int* cols, const double* values) override;
	void ub(int cidx, double val) override;
	void ubs(int cnt, const int* cols, const double* values) override;
	void fixCol(int cidx, double val) override;
	void objcoef(int cidx, double val) override;
	void objcoefs(int cnt, const int* cols, const double* values) override;
	void ctype(int cidx, char val) override;
	void ctypes(int cnt, const int* cols, const char* values) override;
	void switchToLP() override;
	/* Access to underlying XPRESS object */
	XPRSprob getProb() const { return prob; }
private:
	XPRSModel* clone_impl() const override;
	XPRSModel* presolvedmodel_impl() override;
private:
	XPRSprob prob = nullptr;
	bool ownProb = true;
	using SignalHandler = void (*)(int);
	SignalHandler previousHandler = nullptr;
	bool restoreSignalHandler = false;
};

#endif /* XPRSMODEL_H */
