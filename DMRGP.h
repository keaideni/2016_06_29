#ifndef DMRGP_H
#define DMRGP_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "SuperEnergy.h"
#include "Corr.h"
#include <time.h>



class DMRGP
{
public:
	Sub Sys;
	Sub Env;
	Sub newS;
	Sub newE;

	OP truncU;
	OP truncUR;

	Sub m, n;

	QWave fwave;//to store the final wavefunction.
//===============to store the transform wave=========================
        QWave startwave, startwaveL, startwaveR;//to store the ground wave.
        QWave initwave, initwaveL, initwaveR;//to store the initial wave.
        QWave onewave, onewaveL, onewaveR;//to store the onestep wave.
//====================================================================
        //used to define the order of the one step and two step;
        int calnonestepSM, calnonestepSN, calnonestepEM, calnonestepEN;
        int calntwostepSM, calntwostepSN, calntwostepEM, calntwostepEN;
        int caln;





	int OrbitalM, OrbitalN; //the label for the way of growth on blocks. And also the label for the sites M and N.

	int Gdir;//Growth direction: for the absobe way of the two blocks. 1 for the m and -1 for the n.

//=================to caculate the correlation function of the system.===============================
        Corr corr;
        Corr corrn, corrdag;
        double correlation;
        void CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag);
        void CorrUpdate(const int& dir, const Parameter& para);




	double LEnergy, Energy, REnergy;
	OP DenOPWave;

	double FEnergy;
	double FTrace;
	double FTruncerr;
        double FEntanglement;


	double MiuP, MiuN;

	std::ofstream SaveAll;
	std::ofstream Fdata;

	clock_t begin;

	double saveT;


	DMRGP();
	~DMRGP();
	DMRGP(Parameter& para);






	//=================periodic condition==================
	void BuildUpP(Parameter& para, int& OS, int& OE, int& dir);
	void befortruncateP(const Parameter& para, const int& OS, const int& OE) const;
	void getEnergyP(Parameter& para, int dir);
	void truncUpdateP(const Parameter& para, int& OS, int& OE, int dir);


	//================sweep======================
	void SweepP(Parameter& para, int& OS, int& OE, int& dir);
	void getEnergySweepP(Parameter& para, int dir);//dir =1 for the right direction, dir =-1 for the left direction.
	void truncUpdateSweepP(const Parameter& para, int& OS, int& OE, int dir);






};
#endif