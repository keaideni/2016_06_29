#ifndef QWAVE_H
#define QWAVE_H
#include <string>
#include <cmath>
#include <iostream>

#include "Sub.h"

struct classcomm
{
        bool operator()(const std::pair<int, int>& l, const std::pair<int, int>& r) const
        {
                return (l.first + l.second < r.first + r.second) 
                        || (l.first + l.second == r.first + r.second && l.first < r.first);
        }
};




class QWave
{
public:
	//                             m         n        sys is L and env is R
	std::map<std::pair<int, int>, OP, classcomm> WavePart;

	





	QWave();
	~QWave();
	QWave(const QWave& wave);
	void initial(const QWave& wave); //for the inital wave.
	QWave(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot);
	void Initial(const OP& Sys, const OP& m, const OP& n, const OP& Env, int QTot);


	void setZero();
	void normalize();
	int getDim() const;          //calculate the dimention of the wave function.
	void add(const QWave& a, const QWave& b);


	void Wave2f(std::vector<double>& f) const;       //transform the WavePart to a vector, which is used for the calculation of the eigenstate.
	void f2Wave(const std::vector<double>& f);      /*transform the vector to a WavePart, which is used for the calculation of the operator operator on the wave function.*/
	void f2Wave(const VectorXd& f);


	//translate the QWave to OP
	void Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env) const;
	//=======the m in the middle and the n in the right side of the model============
	void Wave2OP(OP& O, const OP& sys, const OP& m, const OP& n, const OP& env, const int& way) const;


	void operator=(const QWave& wave);
	QWave operator+(const QWave& wave);



	void show() const;
	void clear();


	//=============the partition of OPWave======================
	//|WavePart> = O |wave>
	void OPWave2New(const QWave& wave, const OP& O, int flag);    //1 for System, 2 for point M, 3 for point N, 4 for Environment.
	void OSWave2New(const OP& O, const QWave& wave);
	void OEWave2New(const OP& O, const QWave& wave);
	void OMWave2New(const OP& O, const QWave& wave);
	void ONWave2New(const OP& O, const QWave& wave);


	//|storewave> = O|WavePart> + |storewave>==============why it has the second term +|storewave>?==================
	//this part is instructed because we can use it to make less temporary created variation.
	void OPWave(QWave& storewave, const OP& O, int flag) const;     //1 for System, 2 for point M, 3 for point N, 4 for Environment.
	void OSWave(const OP& O, QWave& storewave) const;
	void OEWave(const OP& O, QWave& storewave) const;
	void OMWave(const OP& O, QWave& storewave) const;
	void ONWave(const OP& O, QWave& storewave) const;



	//==============the operation operate on the QWave============================!!!!!!!!!! we can use the reload function to rewrite here.
	//|storewave> = OPl OPr |wavePart> +|storewave>
	/*void BlockWave(int l, const OP& OPl, int r, const OP& OPr, QWave& storewave) const;

	//|storewave>=OPl OPm OPr|PartWave> + |storewave>
	void BlockWave(int l, const OP& OPl, int m, const OP& OPm, int r, const OP& OPr, QWave& storewave) const;

	//|storewave>=OPl OPm OPn OPr|PartWave> + |storewave>
	void BlockWave(int l, const OP& OPl, int m, const OP& OPm, int n, const OP& OPn, int r, const OP& OPr, QWave& storewave) const;*/


	//=============for the initialation

	void onestepSM(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSM, const OP& truncEN);
	void onestepSN(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSN, const OP& truncEM);
	void onestepEM(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSM, const OP& truncEN);
	void onestepEN(const QWave& wave, const OP&sys, const OP&m, const OP&Env, const OP& n, const OP& truncSM, const OP& truncEN);
	void twostepSM(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n);
	void twostepSN(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n);
	void twostepEM(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n);
	void twostepEN(const QWave&wave, const OP&Sys, const OP&m, const OP& Env, const OP&n);



};





#endif