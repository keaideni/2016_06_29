#include "DMRGP.h"
#include "Sub.h"
#include "QWave.h"
#include "Super.h"


DMRGP::DMRGP()
{

}
DMRGP::~DMRGP()
{

}

DMRGP::DMRGP(Parameter& para)
{
	clock_t Abegin;
	double allT;
	int OS(1);
	int OE;
	int dir(1);


	//===========this label for the way of growth the blocks======
	OrbitalM = 2;
	OrbitalN = 2;
	Gdir = 1;
	//====================================================



	OE = para.LatticeSize;


	SaveAll.open("./result/SaveAll");
	if (!SaveAll.is_open())
	{
		//std::cout << "the file doesn't exit!" << std::endl;
	}

	saveT = 0;

	SaveAll << "==============build up================" << std::endl;
	Abegin = clock();





	BuildUpP(para, OS, OE, dir);

	allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
	SaveAll << "===========build up finished============" << std::endl;

	SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
	SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;

	SaveAll << "===============Sweep================" << std::endl;
	saveT = 0;
	Abegin = clock();


        Fdata.open("./result/data");
        if (!Fdata.is_open())
        {
                std::cout << "the file doesn't exit!" << std::endl;
        }

        calnonestepSM=calnonestepSN=calnonestepEM=calnonestepEN=0;
        calntwostepSM=calntwostepSN=calntwostepEM=calntwostepEN=1;
        caln = 0;
	SweepP(para, OS, OE, dir);
//=========================calculate the correlation==================================
        std::ofstream outfile1, outfile2;
        if(OrbitalM % 2 == 1)
        {
                                        
                outfile1.open("./result/resonator.txt");
                
                outfile2.open("./result/qubit.txt");
                }else
                {
                                        
                        outfile1.open("./result/resonator.txt");
                                        
                        outfile2.open("./result/qubit.txt");
                }

                int i(OrbitalM + 1);//this label for Sigma;
                int j(OrbitalM - 1);//this label for Sigmadag and N;
                int fflag(1);
                while(true)
                {
                                        
                                        

                        corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                                        //corrn.show();Sys.SubSysEye.show();

                        CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO);

                        int distence(i - j);

                        std::cout<<"Distence = " << distence << ", the correlation = "
                        <<correlation<<std::endl;
                        outfile1<<"Distence = " << distence << ", the correlation = "
                        <<correlation<<std::endl;


                        if(fflag == 1)
                        {
                                i += 2;
                        }else
                        {
                                j-= 2;
                        }
                        if((i>para.LatticeSize/2)||(j<0))break;
                        fflag *= -1;


                }



                i=(OrbitalM + 2);//this label for Sigma;
                j=(OrbitalM - 2);//this label for Sigmadag and N;
                fflag=1;
                while(true)
                {
                                        
                                        

                        corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                        //corrn.show();Sys.SubSysEye.show();

                        CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO);

                        int distence(i - j);

                        std::cout<<"Distence = " << distence << ", the correlation = "
                        <<correlation<<std::endl;
                        outfile2<<"Distence = " << distence << ", the correlation = "
                        <<correlation<<std::endl;


                        if(fflag == 1)
                        {
                                i += 2;
                        }else
                        {
                                j-= 2;
                        }
                        if((i>para.LatticeSize/2)||(j<0))break;
                        fflag *= -1;


                }

                outfile1.close();
                outfile2.close();
//======================================================================================================
	allT = difftime(clock(), Abegin) / CLOCKS_PER_SEC;
	SaveAll << "===========Sweep finished==============" << std::endl;
	SaveAll << "The build up process takes " << allT << " seconds!" << std::endl;
	SaveAll << "The eigenstate calculation takes " << saveT << " seconds!" << std::endl;
	SaveAll.close();

	


	

	Fdata.close();
}










void DMRGP::BuildUpP(Parameter& para, int& OS, int& OE, int& dir)
{

	befortruncateP(para, OS, OE);
	//int judge = (OE-OS ==1)? 1:0;
	int i(1);

	while (true)
	{
		SaveAll << i << std::endl;

		Sys.read(OS);//Sys.show();
		Env.read(OE);//Env.show();

		m.Initial(para, OrbitalM);
		n.Initial(para, OrbitalN);

		getEnergyP(para, dir);






		if ((OE - OS) == 3)break;


		truncUpdateP(para, OS, OE, dir);

		if (Gdir == 1)
		{
			OrbitalM += 1;
		}
		else
		{
			OrbitalN += 1;
		}
		Gdir *= -1;


		i++;
	}

}


void DMRGP::befortruncateP(const Parameter& para, const int& OS, const int& OE) const
{
	Sub Sys1(para, OS);
	Sys1.save();
	//std::cout<<Sys.SubSys.Max<<std::endl;
	Sub Env1(para, 1);
	Env1.Orbital = OE;
	Env1.save();
}



void DMRGP::getEnergyP(Parameter& para, int dir)
{



	int qtot = Sys.Orbital + 1;
	//std::cout<<qtot<<std::endl;
	Super Sup(para, Sys, m, n, Env, qtot);
	//std::cout<<"hehe"<<std::endl;
	begin = clock();
	SuperEnergy Supp(para, Sup);
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
	//std::cout<<"haha"<<std::endl;
	double trace;
	double truncerr;


	OP temp;
	OP dentemp;


	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);

	dentemp.getDenS(temp);
	truncU.DengetTruncU(para, dentemp, trace, truncerr);//temp.show();

	temp.clear();
	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, -1 * Gdir);




	//Sup.Wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	dentemp.getDenE(temp);
	truncUR.DengetTruncU(para, dentemp, trace, truncerr);

	//truncU.show();
	SaveAll << "Q = " << qtot << "    WaveD = " << std::setw(4) << Sup.Dim
		<< "      OS =" << std::setw(2) << Sys.Orbital << ",  OE =" << std::setw(2) << Env.Orbital
		<< ",    E = " << std::setw(10) << std::setprecision(15) << para.Energy << ",    trace = " << std::setprecision(15) << trace
		<< ",    truncerr = " << std::setprecision(15) << truncerr << std::endl;


	std::cout << "Q = " << qtot << "    WaveD = " <<std::setw(4)<< Sup.Dim
	<< "      OS ="  <<std::setw(2)<<Sys.Orbital << ",  OE =" <<std::setw(2)<< Env.Orbital
	<< ",    E = " <<std::setw(10)<< std::setprecision(15)<<para.Energy <<",    trace = "<< std::setprecision(15)<<trace
	<<",    truncerr = "<< std::setprecision(15)<<truncerr<<std::endl;
	//FEnergy = para.Energy;
	//FTrace = trace;
	//FTruncerr = truncerr;



}


void DMRGP::truncUpdateP(const Parameter& para, int& OS, int& OE, int dir)
{

	OS += dir;
	OE -= dir;

	if (Gdir == 1)
	{
		if (m.QorRl == 0)
		{
			newS.update(para, OS, Sys, m, para.gr);
			newE.update(para, OE, Env, m, para.gl);

		}
		else
		{
			newS.update(para, OS, Sys, m, para.gl);
			newE.update(para, OE, Env, m, para.gr);
		}

	}
	else
	{
		if (n.QorRr == 0)
		{
			newS.update(para, OS, n, Sys, para.gl);
			newE.update(para, OE, n, Env, para.gr);
		}
		else
		{
			newS.update(para, OS, n, Sys, para.gr);
			newE.update(para, OE, n, Env, para.gl);
		}

	}





	//if(OS != 1)
	//{
	newS.trunc(truncU);
	newE.trunc(truncUR);
	//}

	newS.save();





	newE.save();

        truncU.truncsave(OS);
        truncUR.truncsave(OE);

}





//=============sweep================
void DMRGP::SweepP(Parameter& para, int& OS, int& OE, int& dir)
{

	
	para.read();
	FEnergy = 10000000000;
        for(int i = 0; i < 8; ++i)
        {
                para.D += 50;
                int flag(0);
	while (flag<para.SweepNo)
	{

		SaveAll << "the " << (flag + 1) << "th Sweep" << std::endl;
		std::cout<<"the "<<(flag+1)<<"th Sweep"<<std::endl;
		//dir*=(-1);//local here for the first left direction sweep


		//FEnergy = 1000000000;
		while (true)
		{
			m.Initial(para, OrbitalM);
			n.Initial(para, OrbitalN);
			//===============consider the ways in the xishoudian d fangshi ==================================================================================
			//==========saomiao guocheng zhong d zhangdian, zai zhangdao zhongdian d shihou you yige zhangdian fangxiang d fanzhuang, yuanben yinggai ============
			//==========youbian zhangdian d shihou huancheng l zuobian zhangdian.=========================================================================
			//=======zhuyao kan Wave2OP d fangshi====================


                        Sys.read(OS);//Sys.show();
                        Env.read(OE);//Env.show();

//====================================two step of the initial wave====================================================================
                        if(dir == 1)
                        {
                                if(Gdir == -1)
                                {
                                        if(calnonestepSM == calntwostepSM)
                                        {
                                                //ffwave.show();
                                                initwave.twostepSM(onewave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveL.twostepSM(onewaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveR.twostepSM(onewaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);

                                                //fwave11.show();
                                                //exit(true);
                                                ++calntwostepSM;
                                        }
                                }else
                                {
                                        if(calnonestepSN == calntwostepSN)
                                        {
                                                initwave.twostepSN(onewave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveL.twostepSN(onewaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveR.twostepSN(onewaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);

                                                //fwave11.show();exit(true);
                                                ++calntwostepSN;
                                        }
                                }

                        }else
                        {
                                if(Gdir == 1)
                                {
                                        if(calnonestepEM == calntwostepEM)
                                        {
                                                initwave.twostepEM(onewave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveL.twostepEM(onewaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveR.twostepEM(onewaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);

                                                //fwave11.show();exit(true);
                                                ++calntwostepEM;
                                        }
                                }else
                                {
                                        if(calnonestepEN == calntwostepEN)
                                        {
                                                //ffwave.show();
                                                initwave.twostepEN(onewave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveL.twostepEN(onewaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);
                                                initwaveR.twostepEN(onewaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye);

                                                //fwave11.show();exit(true);
                                                ++calntwostepEN;
                                        }

                                }
                        }
//==================================================this position is very important================================================


			//==================this aprt is for the first right sweep, if first left sweep, it should absent==========
			if (OS == (para.LatticeSize - 2) / 2)
			{
				Gdir *= -1;

			}
			//============================present with line dir *= -1 at the end of while(true)=================================

			

			


			OP truncU;

			//exit(1);

			getEnergySweepP(para, dir);


			


			//this one is for the break point at the middle of the line.

			if((flag == para.SweepNo -1) &&(OS == (para.LatticeSize - 2) / 2))
			{
                                




				flag = para.SweepNo;
				break;
			}

			if (dir == 1)
			{

				if (para.LatticeSize - OS == 3)
				{

					break;
				}
			}
			else
			{
				if (OE == 4)
				{

					break;
				}
			}

			
			truncUpdateSweepP(para, OS, OE, dir);




			if (dir == 1)
			{
				if (Gdir == 1)
				{
					OrbitalM += 1;
				}
				else
				{
					OrbitalN += 1;
				}
				Gdir *= -1;
			}
			else
			{

				if (Gdir == 1)
				{
					OrbitalN -= 1;
				}
				else
				{
					OrbitalM -= 1;
				}
				Gdir *= -1;
			}


			//==================this aprt is for the first left sweep, if first left sweep, it should absent==========
			/*if(OS == (para.LatticeSize -2)/2)
			Gdir *= -1;*/
			//============================present with line dir *= -1 before while(true)========================================

                        ++caln;

		}
		dir *= (-1);    //local the for the first right sweep
		flag++;
	}

        Fdata << "Q = " << para.ParticleNo << "    LatticeSize = " << std::setw(4) << para.LatticeSize << ",      gr = " << std::setw(4) << para.gr << ",    gl = " << std::setw(4) << para.gl
                << ",        MiuP = " << std::setprecision(15) << MiuP << ",    MiuN = " << std::setprecision(15) << MiuN << ",    trace = " << std::setprecision(15) << FTrace
                << ",    truncerr = " << std::setprecision(15) << FTruncerr << "              para.D = "<<std::setprecision(15)<<para.D
                <<"          Entanglment = "<<std::setprecision(15)<<FEntanglement<<std::endl;
        std::cout << "Q = " << para.ParticleNo << "    LatticeSize = " << std::setw(4) << para.LatticeSize << ",      gr = " << std::setw(4) << para.gr << ",    gl = " << std::setw(4) << para.gl
                << ",        MiuP = " << std::setprecision(15) << MiuP << ",    MiuN = " << std::setprecision(15) << MiuN << ",    trace = " << std::setprecision(15) << FTrace
                << ",    truncerr = " << std::setprecision(15) << FTruncerr << "              para.D = "<<std::setprecision(15)<<para.D
                <<"          Entanglment = "<<std::setprecision(15)<<FEntanglement<<std::endl;
        }
}






void DMRGP::getEnergySweepP(Parameter& para, int dir)
{


	int qtot = para.ParticleNo;
	double trace;
	double truncerr;
        double Entanglment;




	Super Sup(para, Sys, m, n, Env, qtot);
        //Sup.Wave.show();

        //if(caln != 0)
        //{
                

        //}

        

	begin = clock();
        SuperEnergy Supp;
        if(caln == 0)//(Sys.Orbital == (para.ParticleNo-1)))
        {
	        Supp.init(para, Sup);
        }else
        {

        
                //std::cout<<"Sys.Orbital"<<Sys.Orbital<<std::endl;
                
                Supp.init(para, Sup, initwave);//exit(true); 
                 

        }
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
        startwave = Supp.wave;//temp now



	OP temp;
	OP Dentemp;
	Supp.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	DenOPWave.time(Dentemp, 1.0 / 3);

	Energy = para.Energy;
	para.Energy = 0;
	Super SupL(para, Sys, m, n, Env, qtot - 1);
	begin = clock();
        SuperEnergy SuppL;
        if(caln == 0)//(Sys.Orbital == (para.ParticleNo-1)))
        {
                SuppL.init(para, SupL);
        }else
        {
                SuppL.init(para, SupL, initwaveL);
        }
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
        startwaveL = SuppL.wave;


	temp.clear();
	Dentemp.clear();
	SuppL.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	Dentemp.time(1.0 / 3);
	DenOPWave.addWave(Dentemp);


	LEnergy = para.Energy;


	para.Energy = 0;
	Super SupR(para, Sys, m, n, Env, qtot + 1);
	begin = clock();
        SuperEnergy SuppR;
        if(caln == 0)//(Sys.Orbital == (para.ParticleNo-1)))
        {
	       SuppR.init(para, SupR);
        }else
        {
               SuppR.init(para, SupR, initwaveR);

        }
	saveT += difftime(clock(), begin) / CLOCKS_PER_SEC;
        startwaveR = SuppR.wave;


	temp.clear();
	Dentemp.clear();
	SuppR.wave.Wave2OP(temp, Sys.SubSysEye, m.SubSysEye, n.SubSysEye, Env.SubSysEye, Gdir);
	if (dir == 1)
	{
		Dentemp.getDenS(temp);
	}
	else
	{
		Dentemp.getDenE(temp);
	}
	Dentemp.time(1.0 / 3);
	DenOPWave.addWave(Dentemp);


	REnergy = para.Energy;





	truncU.DengetTruncU(para, DenOPWave, trace, truncerr, Entanglment);//temp.show();






	//truncU.show();
	SaveAll << "Q = " << qtot << ",    E = " << std::setprecision(15) << Energy << ",    LE = " << std::setprecision(15) << LEnergy << ",    RE = " << std::setprecision(15) << REnergy << std::endl
		<< "      OS =" << std::setw(2) << Sys.Orbital << ",  OE =" << std::setw(2) << Env.Orbital
		<< ",    trace = " << std::setprecision(15) << trace
		<< ",    truncerr = " << std::setprecision(15) << truncerr << std::endl << std::endl;




	std::cout << "Q = " << qtot << ",    E = " << std::setprecision(15)<<Energy << ",    LE = " << std::setprecision(15)<<LEnergy  << ",    RE = " << std::setprecision(15)<<REnergy<<std::endl
	<< "      OS ="  <<std::setw(2)<<Sys.Orbital << ",  OE =" <<std::setw(2)<< Env.Orbital
        << "      OrbitalM ="  <<std::setw(2)<<OrbitalM << ",  OrbitalN =" <<std::setw(2)<< OrbitalN
	<<",    trace = "<< std::setprecision(15)<<trace
	<<",    truncerr = "<< std::setprecision(15)<<truncerr << std::endl<<std::endl;

	if (Sys.Orbital == (para.LatticeSize - 2) / 2)
	{
		FEnergy = Energy;
		MiuP = REnergy - Energy;
		MiuN = Energy - LEnergy;
		FTrace = trace;
		FTruncerr = truncerr;
                FEntanglement = Entanglment;
                fwave = Supp.wave;
		
	}



}



void DMRGP::truncUpdateSweepP(const Parameter& para, int& OS, int& OE, int dir)
{

        CorrUpdate(dir, para);
        /*std::cout<<"the corr Dim: "<<std::endl;
        for(auto it = corr.CorrO.QDim.begin(); it != corr.CorrO.QDim.end(); ++it)
        {
                std::cout << it->first << "=>" << it->second<<std::endl;
        }*/

	


	//std::cout<<dir<<std::endl;
	if (dir == 1)
	{

		if (Gdir == 1)
		{
			if (m.QorRl == 0)
			{
				newS.update(para, OS + 1, Sys, m, para.gr);
			}
			else
			{
				newS.update(para, OS + 1, Sys, m, para.gl);
			}

//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncE;
                        truncE.truncread(OE);//truncE.show();
                        onewave.onestepSM(startwave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);
                        onewaveL.onestepSM(startwaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);
                        onewaveR.onestepSM(startwaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);

                        ++calnonestepSM;
                        //fwave1.show();ffwave.show();exit(true);
//============================================================================================================================================================
		}
		else
		{
			if (n.QorRr == 0)
			{
				newS.update(para, OS + 1, n, Sys, para.gl);
			}
			else
			{
				newS.update(para, OS + 1, n, Sys, para.gr);
			}
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncE;
                        truncE.truncread(OE);//truncE.show();
                        onewave.onestepSN(startwave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);
                        onewaveL.onestepSN(startwaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);
                        onewaveR.onestepSN(startwaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncE);

                        ++calnonestepSN;
                        //ffwave.show();//exit(true);
//=============================================================================================================================================================
		}
		//newS.show();
		//truncU.show();

		newS.trunc(truncU);
		//newS.show();

                truncU.truncsave(OS+1);

		newS.save();
	}
	else
	{
		


		if (Gdir == 1)
		{
			if (n.QorRr == 0)
			{
				newS.update(para, OE - 1, n, Env, para.gr);
			}
			else
			{
				newS.update(para, OE - 1, n, Env, para.gl);
			}
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncS;
                        truncS.truncread(OS);//truncE.show();
                        onewave.onestepEN(startwave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);
                        onewaveL.onestepEN(startwaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);
                        onewaveR.onestepEN(startwaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);

                        ++calnonestepEN;
                        //ffwave.show();exit(true);
//============================================================================================================================================================
		}
		else
		{


			if (m.QorRl == 0)
			{
				newS.update(para, OE - 1, Env, m, para.gl);
			}
			else
			{
				newS.update(para, OE - 1, Env, m, para.gr);
			}
//============================================================================================================================================================
                        //this part is for the wavetransform.
                        //QWave ffwave;
                        OP truncS;
                        truncS.truncread(OS);//truncE.show();
                        onewave.onestepEM(startwave, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);
                        onewaveL.onestepEM(startwaveL, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);
                        onewaveR.onestepEM(startwaveR, Sys.SubSysEye, m.SubSysEye, Env.SubSysEye, n.SubSysEye, truncU, truncS);

                        ++calnonestepEM;
                        //ffwave.show();exit(true);
//============================================================================================================================================================
			
		}




		newS.trunc(truncU);
		//truncUR.show();
                truncU.truncsave(OE-1);

		newS.save();
	}

        /*std::cout<<"the newS Dim: "<<std::endl;
        for(auto it = newS.SubSysEye.QDim.begin(); it != newS.SubSysEye.QDim.end(); ++it)
        {
                std::cout << it->first << "=>" << it->second<<std::endl;
        }*/
	OE += dir;
	OS += dir;

        //fwave.show();
        
}




void DMRGP::CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag)
{
	QWave wave1, wave2;
	std::vector<double> f, f1;
	fwave.Wave2f(f);
        QWave ffwave(fwave);
	double number(0);
	wave1.clear();
	wave1.OSWave2New(corrn, fwave);
        ffwave.initial(wave1);
	ffwave.Wave2f(f1);

	for(int i = 0; i < f.size(); ++i)
	{
		number += f[i]*f1[i];
	}
	wave1.clear();
	wave1.OSWave2New(corrcdag, fwave);
	wave2.clear();
	wave2.OEWave2New(corrc, wave1);
        ffwave.initial(wave2);
	ffwave.Wave2f(f1);
	correlation = 0;
	for(int i = 0; i < f.size(); ++i)
	{
		correlation += f[i]*f1[i];
	}
	correlation /= number;


}

void DMRGP::CorrUpdate(const int& dir, const Parameter& para)
{
        //this is an interesting phenomenon, why this part must be in this position?!
        if((dir == -1)&&(OrbitalM > para.LatticeSize/4)&&(OrbitalN > para.LatticeSize/4 + 1 ) )
        {
                //when the block Env eat the point m, the corr1 begins to initialize and update.
                if(OrbitalM == OrbitalN)
                {
                                        
                        //if(OrbitalN == para.LatticeSize/4 + 1) return;
                        corr.Initial(Env, m, OrbitalM, 1, truncU);
                        corr.save();
                        for(int i = para.LatticeSize/2; i > OrbitalM; --i)
                        {
                                corr.read(i, 1);
                                corr.update(m, truncU, 1);
                                corr.save();
                        }
                }//when the block Env eat the point n, the corr1 only update.
                else
                {
                        for(int i = para.LatticeSize/2; i >= OrbitalN; --i)
                        {
                                corr.read(i, 1);
                                corr.update(n, truncU, 2);
                                corr.save();
                        }
                }

        }else if((dir == 1) && (OrbitalM <= (para.LatticeSize ) / 4+1 )&&(OrbitalN<=para.LatticeSize/4))
        {
                //when the block Sys eat the point m, the corr2 and corr3 begin to initialize and update.
                if(OrbitalM == OrbitalN)
                {
                        corr.Initial(Sys, m, OrbitalM, 2, truncU);
                        corr.save();
                        corr.Initial(Sys, m, OrbitalM, 3, truncU);
                        corr.save();
                        for(int i = 2; i < OrbitalM; ++i)
                        {
                                corr.read(i, 2);
                                corr.update(m, truncU, 1);
                                corr.save();


                                corr.read(i, 3);
                                corr.update(m, truncU, 1);
                                corr.save();

                        }
                }//when the block Env eat the point n, the corr1 only update.
                else
                {
                        for(int i = 2; i <= OrbitalN; ++i)
                        {
                                corr.read(i, 2);
                                //corr.corro().show();
                                /*OP temp;
                                temp.kronO(n.SubSysEye, corr.corro());*/
                                corr.update(n, truncU, 2);
                                corr.save();

                                corr.read(i, 3);
                                
                                //corr.corro().show();
                                
                                corr.update(n, truncU, 2);
                                corr.save();


                        }
                }
        }
}