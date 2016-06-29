#include "Super.h"
#include <ctime>


Super::Super(){}


Super::~Super()
{

}






Super::Super(const Parameter& para_, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, const int& TotQ)
{
	Wave.Initial(sys.SubSys, m.SubSys, n.SubSys, env.SubSys, TotQ);
	SaveWave = Wave;

	para = para_;

	Dim = Wave.getDim();

	RealSysBlock = &sys;
	RealMBlock = &m;
	RealNBlock = &n;
	RealEnvBlock = &env;
}



void Super::Initial(const Parameter& para_, const Sub& sys, const Sub& m, const Sub& n, const Sub& env, int TotQ)
{
	Wave = QWave(sys.SubSys, m.SubSys, n.SubSys, env.SubSys, TotQ);
	SaveWave = Wave;

	para = para_;

	Dim = Wave.getDim();

	RealSysBlock = &sys;
	RealMBlock = &m;
	RealNBlock = &n;
	RealEnvBlock = &env;
}






//the point f translate to g
void Super::f1tof2(const double* f, double* g)
{

	std::vector<double> f1, g1;
	for (int i = 0; i<Dim; i++)
	{
		f1.push_back(f[i]);
	}

	f1tof2(f1, g1);

	for (int i = 0; i<Dim; i++)
	{
		g[i] = g1[i];
	}
}



//the vector f translate to g.
void Super::f1tof2(const std::vector<double>& f, std::vector<double>& g)
{

	g.clear();
	Wave.f2Wave(f);
	SaveWave.setZero();
	OneStep(); //for the |SaveWave> = H |Wave>
	SaveWave.Wave2f(g);
}






//===========to calculate the QWave after the operation of the hamiltion=============
void Super::OneStep()
{

	OP p1, p2;
	//=========the term have no relationship to the boundary condition========
	Wave.OPWave(SaveWave, RealSysBlock->SubSys, 1);//sys

	Wave.OPWave(SaveWave, RealMBlock->SubSys, 2);//M

	Wave.OPWave(SaveWave, RealNBlock->SubSys, 3);//N

	Wave.OPWave(SaveWave, RealEnvBlock->SubSys, 4);//env




	if (RealMBlock->QorRl == 0)
	{
		//Sys-M
		p1.time(para.gr, RealSysBlock->SubSysC);
		p2.time(para.gr, RealSysBlock->SubSysCdag);
		BlockWave(1, p1, 2, RealMBlock->SubSysCdag);
		BlockWave(1, p2, 2, RealMBlock->SubSysC);



		//========the term depend on the boundary condition==============
		if (para.EdgeCondition == 0)//open condition
		{
			//M-N
			p1.time(para.gl, RealMBlock->SubSysC);
			p2.time(para.gl, RealMBlock->SubSysCdag);

			BlockWave(2, p1, 3, RealNBlock->SubSysCdag);
			BlockWave(2, p2, 3, RealNBlock->SubSysC);

			//N-Env
			p1.time(para.gr, RealNBlock->SubSysC);
			p2.time(para.gr, RealNBlock->SubSysCdag);

			BlockWave(3, p1, 4, RealEnvBlock->SubSysCdag);
			BlockWave(3, p2, 4, RealEnvBlock->SubSysC);


		}
		else//==================periodic condition=====================
		{
			//M-Env
			p1.time(para.gl, RealMBlock->SubSysC);
			p2.time(para.gl, RealMBlock->SubSysCdag);

			BlockWave(2, p1, 4, RealEnvBlock->SubSysCdag);
			BlockWave(2, p2, 4, RealEnvBlock->SubSysC);


			if (RealNBlock->QorRl == 0)
			{

				//Env-N
				p1.time(para.gr, RealEnvBlock->SubSysC1);
				p2.time(para.gr, RealEnvBlock->SubSysCdag1);

				BlockWave(4, p1, 3, RealNBlock->SubSysCdag);
				BlockWave(4, p2, 3, RealNBlock->SubSysC);

				//N-Sys
				p1.time(para.gl, RealNBlock->SubSysC);
				p2.time(para.gl, RealNBlock->SubSysCdag);

				BlockWave(3, p1, 1, RealSysBlock->SubSysCdag1);
				BlockWave(3, p2, 1, RealSysBlock->SubSysC1);
			}
			else
			{

				//Env-N
				p1.time(para.gl, RealEnvBlock->SubSysC1);
				p2.time(para.gl, RealEnvBlock->SubSysCdag1);

				BlockWave(3, RealNBlock->SubSysCdag, 4, p1);
				BlockWave(3, RealNBlock->SubSysC, 4, p2);

				//N-Sys
				p1.time(para.gr, RealNBlock->SubSysC);
				p2.time(para.gr, RealNBlock->SubSysCdag);

				BlockWave(3, p1, 1, RealSysBlock->SubSysCdag1);
				BlockWave(3, p2, 1, RealSysBlock->SubSysC1);
			}
		}

	}
	else
	{

		//Sys-M
		p1.time(para.gl, RealSysBlock->SubSysC);
		p2.time(para.gl, RealSysBlock->SubSysCdag);
		BlockWave(1, p1, 2, RealMBlock->SubSysCdag);
		BlockWave(1, p2, 2, RealMBlock->SubSysC);



		//========the term depend on the boundary condition==============
		if (para.EdgeCondition == 0)
		{
			//M-N
			p1.time(para.gr, RealMBlock->SubSysC);
			p2.time(para.gr, RealMBlock->SubSysCdag);

			BlockWave(2, p1, 3, RealNBlock->SubSysCdag);
			BlockWave(2, p2, 3, RealNBlock->SubSysC);

			//N-Env
			p1.time(para.gl, RealNBlock->SubSysC);
			p2.time(para.gl, RealNBlock->SubSysCdag);

			BlockWave(3, p1, 4, RealEnvBlock->SubSysCdag);
			BlockWave(3, p2, 4, RealEnvBlock->SubSysC);
		}
		else                                                  //for periodic condition
		{
			//M-Env
			p1.time(para.gr, RealMBlock->SubSysC);
			p2.time(para.gr, RealMBlock->SubSysCdag);

			BlockWave(2, p1, 4, RealEnvBlock->SubSysCdag);
			BlockWave(2, p2, 4, RealEnvBlock->SubSysC);

			if (RealNBlock->QorRl == 0)
			{


				//Env-N
				p1.time(para.gr, RealEnvBlock->SubSysC1);
				p2.time(para.gr, RealEnvBlock->SubSysCdag1);

				BlockWave(4, p1, 3, RealNBlock->SubSysCdag);
				BlockWave(4, p2, 3, RealNBlock->SubSysC);

				//N-Sys
				p1.time(para.gl, RealNBlock->SubSysC);
				p2.time(para.gl, RealNBlock->SubSysCdag);

				BlockWave(3, p1, 1, RealSysBlock->SubSysCdag1);
				BlockWave(3, p2, 1, RealSysBlock->SubSysC1);
			}
			else
			{

				//Env-N
				p1.time(para.gl, RealEnvBlock->SubSysC1);
				p2.time(para.gl, RealEnvBlock->SubSysCdag1);

				BlockWave(4, p1, 3, RealNBlock->SubSysCdag);
				BlockWave(4, p2, 3, RealNBlock->SubSysC);

				//N-Sys
				p1.time(para.gr, RealNBlock->SubSysC);
				p2.time(para.gr, RealNBlock->SubSysCdag);

				BlockWave(3, p1, 1, RealSysBlock->SubSysCdag1);
				BlockWave(3, p2, 1, RealSysBlock->SubSysC1);
			}
		}
	}


	

}



void Super::BlockWave(const int& l, const OP& OPl, const int& r, const OP& OPr)
{
	tempWave.clear();

	tempWave.OPWave2New(Wave, OPr, r);
	tempWave.OPWave(SaveWave, OPl, l);

}



//============translate the base state to Qwave form==============
void Super::normalizedCopy(const double* f0)
{
	std::vector<double> v;
	double norm(0);

	for (int i = 0; i<Dim; i++)
	{
		norm += f0[i] * f0[i];
	}

	norm = sqrt(norm);
	for (int i = 0; i<Dim; i++)
	{
		v.push_back(f0[i] / norm);
	}

	Wave.f2Wave(v);

}




void Super::show()
{
	std::cout << "the System part: " << std::endl;
	RealSysBlock->show();
	std::cout << "the M part: " << std::endl;
	RealMBlock->show();
	std::cout << "the N part: " << std::endl;
	RealNBlock->show();
	std::cout << "the Env part: " << std::endl;
	RealEnvBlock->show();
	std::cout << "the QWave part: " << std::endl;
	Wave.show();
	std::cout << "Dim = " << Dim << std::endl;
}