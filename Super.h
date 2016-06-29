#ifndef SUPER_H
#define SUPER_H

#include <iostream>
#include <vector>
#include <string>
#include "QWave.h"

class Super
{
private:
	double sigma_;
public:


	int rows() { return Dim; };
	int cols() { return Dim; };
	void set_shift(double sigma) { sigma_ = sigma; }
	void perform_op(double *x_in, double *y_out)
	{
		f1tof2(x_in, y_out);
	};
	



	const Sub* RealSysBlock;
	const Sub* RealEnvBlock;
	const Sub* RealMBlock;
	const Sub* RealNBlock;


	QWave Wave, tempWave;
	QWave SaveWave;
	Parameter para;
	int Dim;

	



	Super();
	~Super();
	Super(const Parameter& para, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, const int& TotQ);

	void Initial(const Parameter& para, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, int TotQ);

	//===========to calculate the QWave after the operation of the hamiltion=============
	void f1tof2(const double* f, double* g);
	void f1tof2(const std::vector<double>& f, std::vector<double>& g);
	void OneStep();


	void BlockWave(const int& l, const OP& OPl, const int& r, const OP& OPr);


	void normalizedCopy(const double* f0);




	void show();

};











#endif