#ifndef PARAMETER_H
#define PARAMETER_H
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
class Parameter
{

public:

	double Wz;
	double Wc;
	double gr;
	double gl;

	double Energy;

	int LatticeSize;
	int ParticleNo;     //total number of the whole system;
	int SiteNo;        //the number set for a ceter, and some number fluent aroud it.
	int DeltaQL;       //the number in the left of the ceter number.
	int DeltaQR;       //the nmmber in the right of the ceter number.

	int D;
	int SweepNo;
	int EdgeCondition; //0 for open boundary condition.




	Parameter();
	~Parameter();
	Parameter(const Parameter& para);



	void save();
	void read();
	void show();

};






















#endif
