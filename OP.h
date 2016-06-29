#ifndef OP_H
#define OP_H

#include <string>
#include <utility>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Eigenvalues> 
#include <Eigen/Dense>
#include "Parameter.h"



struct classcom
{
	size_t operator()(const std::pair<int, int>& l) const
	{
		return (size_t) ((l.first<<20) ^l.second);
	}
};


using namespace Eigen;

class OP
{
public:
	std::unordered_map<int, int> QDim;
	std::unordered_map<int, MatrixXd> QMat;
	std::unordered_map<int, int> RLQ;


	static int Max;

	OP();
	~OP();
	OP(const int& Dmin, const int& Dmax, const int& Qmin, const int& Qmax, const int& delta, const int& str);

	OP(const OP& a);


	void iniRaiseQ(const int& Dmin, const int& Dmax, const int& Qmin, const int& Qmax, const int& delta, const int& str);
	//1.for boson annihilation. 2.for boson creation. 3 for I. 4 for N.
	void iniRaiseQ(const int str);        
	//used for creat a qubit operator: case 1 for sigmaZ, case 2 for sigmadag, case 3 for sigma, case 4 for Eye.
	void getValue(const int& str, const int& i, double& y);


	void findDim(const OP& a, const OP& b, std::unordered_map<int, int>& oldDim, std::unordered_map<std::pair<int, int>, int, classcom>& startDim) const;
	void kronO(const OP& a, const OP&b);
	//void getmat(MatrixXd& oldmat, const MatrixXd& tempmat, const int& startL1, const int& startR1);
	void transO(const OP& a);


	void add(const OP& a, const OP& b);               //for the operator +, but don't need to copy anything.
	void add(const OP& a);                 //fot the operator +=, but don't need to copy anything.
	void time(const double& x);
	void time(const OP& a, const double& x);
	void time(const double& x, const OP& a);
	void time(const OP& a, const OP& b);


	//for the QWave OP
	int ltime(const OP& a, const OP& wave);      
	//this function is used for the operator a function on the wave operator wave. 
	//the return of this function is used to see whether the time is finished or not.
	int rtime(const OP& a, const OP& wave);      
	//this two functions are used for the operator of System(the a in l), and Envirment(the a in r)

	//wave operator, for which don't have the QDim. (the operator in the QWave can't have a QDim)
	void addWave(const OP& a, const OP& b);
	void addWave(const OP& a);


	//for the truncation
	void getTruncU(const Parameter& para, const OP& OPWave);
	void getTruncUR(const Parameter& para, const OP& OPWave);
	void getTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
	void getTruncUR(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
	void truncL(const OP& trunc, const OP& O);
	void truncR(const OP& O, const OP& trunc);
	void trunc(const OP& truncU);


	void getDenS(const OP& OPWave);
	void getDenE(const OP& OPWave);
	void DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
        void DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr, double& Entanglement);



	OP operator+(const OP& a);
	OP operator*(const double& x);
	OP operator*(const OP& a);
	void operator=(const OP& a);




	void show() const;
	void clear();
	void save(std::ofstream& outfile);
	void read(std::ifstream& infile);



        void truncsave(const int& orbital);
        void truncread(const int& orbital);
};



OP operator*(const double& x, OP& a);




std::string itos(int i);














#endif