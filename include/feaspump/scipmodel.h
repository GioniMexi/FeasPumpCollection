/**
 * @file scipmodel.h
 * @brief Implementation of MIPModelI for SCIP
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
 */

#ifndef SCIPMODEL_H
#define SCIPMODEL_H

#include "mipmodel.h"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/scip_copy.h>

class SCIPModel: public MIPModelI
{


public:
	SCIPModel();
	SCIPModel(SCIP* scip, bool _ownProb = false);
	~SCIPModel() override;
	std::unique_ptr<SCIPModel> clone() const { return std::unique_ptr<SCIPModel>(this->clone_impl()); }
	/* Read/Write */
	void readModel(const std::string& filename) override;
	void writeModel(const std::string& filename, const std::string& format="") const override;
	void writeSol(const std::string& filename) const override;
	/* Solve */
	double lpopt(char method, bool decrease_tol, bool initial) override;
	void mipopt() override;
	/* Presolve/postsolve */
	void presolve() override;
	void postsolve() override;
	std::unique_ptr<SCIPModel> presolvedModel() { return std::unique_ptr<SCIPModel>(this->presolvedmodel_impl()); }
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
	void terminationReason(std::string& reason) override;
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
	void switchToMIP() override;

	/* Access to underlying SCIP object */
	SCIP* getProb() const { return scip; }

	std::vector<std::string> getVarNames();

private:
	SCIPModel* clone_impl() const override;
	SCIPModel* presolvedmodel_impl() override;
	void copyLpi(SCIP_LPI** targetlpi, SCIP_LPI* srclpi);
	SCIP_VAR ** varsData();
	SCIP_STAGE getProbStage();

private:
	SCIP* scip = nullptr; 
	std::vector<std::string> origColNames;
	std::vector<std::string> transColNames;
	std::vector<std::string> presColNames;
	mutable SCIP_LPI* lpi = nullptr; 
	SCIP_RETCODE retcode;
	SCIP_STAGE stage;
	mutable SCIP_VAR **varsAtStart;

 	bool ownProb = true;
	bool stageLP = false;
	bool isClone = false;
	bool postsolved = false;

	int nvarsStartStage3;
	int nconsStartStage3;

	using SignalHandler = void (*)(int);
	SignalHandler previousHandler = nullptr;
	bool restoreSignalHandler = false;
};


#endif /* SCIPMODEL_H */
