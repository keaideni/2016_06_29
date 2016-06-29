#include "Parameter.h"


Parameter::~Parameter(){}

Parameter::Parameter(){}


Parameter::Parameter(const Parameter &para)
{
	Wz = para.Wz;
	Wc = para.Wc;
	gr = para.gr;
	gl = para.gl;

	Energy = para.Energy;
	LatticeSize = para.LatticeSize;
	ParticleNo = para.ParticleNo;
	SiteNo = para.SiteNo;
	DeltaQL = para.DeltaQL;
	DeltaQR = para.DeltaQR;

	D = para.D;
	SweepNo = para.SweepNo;
	EdgeCondition = para.EdgeCondition;

}

void Parameter::save()
{
	std::ofstream outfile;
	outfile.open("./data/QNosave.txt");
	if (!outfile.is_open())
	{
		std::cout << "the file doesn't exit!" << std::endl;
	}
	outfile << "Wi= " << Wz << std::endl
		<< "Wc= " << Wc << std::endl
		<< "gr= " << gr << std::endl
		<< "gl=" << gl << std::endl
		<< "Energy= " << Energy << std::endl
		<< "LatticeSize= " << LatticeSize << std::endl
		<< "PaticleNo= " << ParticleNo << std::endl
		<< "SiteNo= " << SiteNo << std::endl
		<< "DeltaQL= " << DeltaQL << std::endl
		<< "DeltaQR= " << DeltaQR << std::endl
		<< "D= " << D << std::endl
		<< "SweepNop= " << SweepNo << std::endl
		<< "EdgeCondition= " << EdgeCondition << std::endl;
	outfile.close();
}



void Parameter::read()
{
	std::ifstream infile("./data/QNosave.txt");
	if (!infile.is_open())
	{
		std::cout << "the data file doesn't exit!" << std::endl;
	}
	std::string str;
	infile >> str >> Wz
		>> str >> Wc
		>> str >> gr
		>> str >> gl
		>> str >> Energy
		>> str >> LatticeSize
		>> str >> ParticleNo
		>> str >> SiteNo
		>> str >> DeltaQL
		>> str >> DeltaQR
		>> str >> D
		>> str >> SweepNo
		>> str >> EdgeCondition;
	infile.close();

}



void Parameter::show()
{
	std::cout << "Wi= " << Wz << std::endl
		<< "Wc= " << Wc << std::endl
		<< "gr= " << gr << std::endl
		<< "gl= " << gl << std::endl
		<< "Energy= " << Energy << std::endl
		<< "LatticeSize= " << LatticeSize << std::endl
		<< "PaticleNo= " << ParticleNo << std::endl
		<< "SiteNo= " << SiteNo << std::endl
		<< "DeltaQL= " << DeltaQL << std::endl
		<< "DeltaQR= " << DeltaQR << std::endl
		<< "D= " << D << std::endl
		<< "SweepNop= " << SweepNo << std::endl
		<< "EdgeCondition= " << EdgeCondition << std::endl;
}
