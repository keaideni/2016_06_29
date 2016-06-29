#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <iostream>
#include "Super.h"

using namespace Spectra;

#ifndef SUPERENERGY_H
#define SUPERENERGY_H
class SuperEnergy
{
public:
	QWave wave;

        SuperEnergy(){};
	SuperEnergy(Parameter&para,Super& sup)
	{
		wave = sup.Wave;
		SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, 6);
		eigs.init();
		eigs.compute();
		if (eigs.info() == SUCCESSFUL)
		{
			wave.f2Wave(eigs.eigenvectors(1));
			para.Energy = eigs.eigenvalues()(0);
			std::cout << eigs.num_iterations() << std::endl;
		}

		
	};
        void init(Parameter&para,Super& sup)
        {
                wave = sup.Wave;
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, 6);
                eigs.init();
                eigs.compute();
                if (eigs.info() == SUCCESSFUL)
                {
                        wave.f2Wave(eigs.eigenvectors(1));
                        para.Energy = eigs.eigenvalues()(0);
                        std::cout << eigs.num_iterations() << std::endl;
                }

                
        };

        SuperEnergy(Parameter&para,Super& sup, const QWave& initwave)
        {
                wave = sup.Wave;
                std::vector<double> f;
                wave.initial(initwave);
                wave.Wave2f(f);
                double *pt = new double [sup.Dim];
                for(int i = 0; i < sup.Dim; ++i)pt[i] = f.at(i);
                
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, 6);
                eigs.init(pt);
                eigs.compute();
                if (eigs.info() == SUCCESSFUL)
                {
                        wave.f2Wave(eigs.eigenvectors(1));
                        para.Energy = eigs.eigenvalues()(0);
                        std::cout << eigs.num_iterations() << std::endl;
                }

                
        };
        void init(Parameter&para,Super& sup, const QWave& initwave)
        {
                wave = sup.Wave;
                std::vector<double> f;
                wave.initial(initwave);
                wave.Wave2f(f);
                double *pt = new double [sup.Dim];
                for(int i = 0; i < sup.Dim; ++i)pt[i] = f.at(i);
                
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, 6);
                eigs.init(pt);
                eigs.compute();
                if (eigs.info() == SUCCESSFUL)
                {
                        wave.f2Wave(eigs.eigenvectors(1));
                        para.Energy = eigs.eigenvalues()(0);
                        std::cout << eigs.num_iterations() << std::endl;
                }

                
        };
	

};







#endif
