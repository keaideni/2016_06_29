#ifndef SUB_H_
#define SUB_H_
#include "OP.h"



class Sub
{
public:


	int Orbital;

	//=======these amounts are used for label on two sides of the blocks=========
	int QorRl;//1 for Qubit and 0 for resonator;
	int QorRr;//1 for Qubit and 0 for resonator;
	//==================================================================


	OP SubSys;
	OP SubSysCdag; //2
	OP SubSysC; //1
	OP SubSysEye;  //3



	OP SubSysCdag1;
	OP SubSysC1;



	Sub();
	~Sub();
	Sub(const Sub& block);
	Sub(const Parameter& para, const int& orbital_);
	//=====================================================================================
	//update the whole Sub. the constant coup is for the couple constant.
	Sub(const Parameter& para, const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup);
	//======================================================================================


	void Initial(const Parameter& para, const int & orbital_);
	void update(const Parameter& para, const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup);
	//void update(const Parameter& para, const Sub& old, const int orbital_, const int& way);//for the periodic condition.
	//truncated U to update the Sub.
	void trunc(const OP& truncU);


	void updateBlockH(const Sub& oldL, const Sub& oldR, const double& coup);//update the SubSys
	//void updateBlockH(OP& NewSys, const Parameter& para, const Sub& old, const int& way);//for the periodic condition.



	void operator=(const Sub& block);



	void show() const;
	void clear();
	void save();
	void read(int orbital_);
};





#endif