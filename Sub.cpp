#include "Sub.h"


Sub::Sub(){}

Sub::~Sub(){}


Sub::Sub(const Sub& block)
{
	Orbital = block.Orbital;
	QorRl = block.QorRl;
	QorRr = block.QorRr;
	SubSys = block.SubSys;
	SubSysCdag = block.SubSysCdag;
	SubSysC = block.SubSysC;
	SubSysEye = block.SubSysEye;


	SubSysCdag1 = block.SubSysCdag1;
	SubSysC1 = block.SubSysC1;




}



Sub::Sub(const Parameter& para, const int& orbital_)
{
	Orbital = orbital_;
	QorRl = orbital_ % 2;
	QorRr = orbital_ % 2;

	if (QorRl == 0)
	{
		//std::cout<<"haha"<<std::endl;
		int dmin = para.SiteNo - para.DeltaQL < 0 ? 0 : para.SiteNo - para.DeltaQL;
		int dmax = para.SiteNo + para.DeltaQR > para.ParticleNo ? para.ParticleNo : para.SiteNo + para.DeltaQR;

		SubSysC.iniRaiseQ(dmin, dmax, dmin + 1, dmax, -1, 1);
		SubSysCdag.iniRaiseQ(dmin, dmax, dmin, dmax - 1, 1, 2);
		SubSysEye.iniRaiseQ(dmin, dmax, dmin, dmax, 0, 3);
		SubSys.iniRaiseQ(dmin, dmax, dmin, dmax, 0, 4);


		SubSys.time(para.Wc);



	}
	else
	{
		SubSys.iniRaiseQ(1);
		SubSysCdag.iniRaiseQ(2);
		SubSysC.iniRaiseQ(3);

		SubSysEye.iniRaiseQ(4);
		SubSys.time(para.Wz / 2);


	}


	SubSysC1 = SubSysC;
	SubSysCdag1 = SubSysCdag;




}





void Sub::Initial(const Parameter& para, const int & orbital_)
{

	clear();
	Orbital = orbital_;
	QorRl = orbital_ % 2;
	QorRr = orbital_ % 2;

	if (QorRl == 0)
	{
		//std::cout<<"haha"<<std::endl;
		int dmin = para.SiteNo - para.DeltaQL < 0 ? 0 : para.SiteNo - para.DeltaQL;
		int dmax = para.SiteNo + para.DeltaQR > para.ParticleNo ? para.ParticleNo : para.SiteNo + para.DeltaQR;

		SubSysC.iniRaiseQ(dmin, dmax, dmin + 1, dmax, -1, 1);
		SubSysCdag.iniRaiseQ(dmin, dmax, dmin, dmax - 1, 1, 2);
		SubSysEye.iniRaiseQ(dmin, dmax, dmin, dmax, 0, 3);
		SubSys.iniRaiseQ(dmin, dmax, dmin, dmax, 0, 4);


		SubSys.time(para.Wc);



	}
	else
	{
		SubSys.iniRaiseQ(1);
		SubSysCdag.iniRaiseQ(2);
		SubSysC.iniRaiseQ(3);

		SubSysEye.iniRaiseQ(4);
		SubSys.time(para.Wz / 2);


	}


	SubSysC1 = SubSysC;
	SubSysCdag1 = SubSysCdag;
}








//update the whole Sub.
Sub::Sub(const Parameter& para, const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup)
{
	Orbital = orbital_;
	QorRl = oldL.QorRl;
	QorRr = oldR.QorRr;

	if ((oldL.QorRr + oldR.QorRl) % 2 == 0)
	{
		std::cout << "WARNING!!!!!!" << std::endl;
		std::cout << "the growth of block is wrong!" << std::endl;
		//exit(1); 
	}


	updateBlockH(oldL, oldR, coup);

	SubSysCdag.kronO(oldL.SubSysEye, oldR.SubSysCdag);
	SubSysC.kronO(oldL.SubSysEye, oldR.SubSysC);
	SubSysEye.kronO(oldL.SubSysEye, oldR.SubSysEye);







	if (para.EdgeCondition != 0)
	{

		SubSysC1.kronO(oldL.SubSysC1, oldR.SubSysEye);

		SubSysCdag1.kronO(oldL.SubSysCdag1, oldR.SubSysEye);

	}

	//the rest is not needed in the superblock.


}





//update the SubSys.
void Sub::updateBlockH(const Sub& oldL, const Sub& oldR, const double& coup)
{

	clear();
	OP tempOP1;
	//OP tempOP2;



	SubSys.kronO(oldL.SubSysEye, oldR.SubSys);
	//tempOP2.time(tempOP1, para.Wc);



	//H kron I.

	//tempOP2.clear();
	tempOP1.kronO(oldL.SubSys, oldR.SubSysEye);
	SubSys.add(tempOP1);



	tempOP1.kronO(oldL.SubSysCdag, oldR.SubSysC1);
	//oldL.SubSysCdag.show();
	//oldR.SubSysC1.show();
	tempOP1.time(coup);
	SubSys.add(tempOP1);


	tempOP1.clear();
	tempOP1.kronO(oldL.SubSysC, oldR.SubSysCdag1);
	tempOP1.time(coup);
	SubSys.add(tempOP1);






}







void Sub::update(const Parameter& para, const int& orbital_, const Sub& oldL, const Sub& oldR, const double& coup)
{

	Orbital = orbital_;
	QorRl = oldL.QorRl;
	QorRr = oldR.QorRr;

	if ((oldL.QorRr + oldR.QorRl) % 2 == 0)
	{
		std::cout << "WARNING!!!!!!" << std::endl;
		std::cout << "the growth of block is wrong!" << std::endl;
		exit(1);
	}


	updateBlockH(oldL, oldR, coup);

	SubSysCdag.kronO(oldL.SubSysEye, oldR.SubSysCdag);
	SubSysC.kronO(oldL.SubSysEye, oldR.SubSysC);
	SubSysEye.kronO(oldL.SubSysEye, oldR.SubSysEye);







	if (para.EdgeCondition != 0)
	{

		SubSysC1.kronO(oldL.SubSysC1, oldR.SubSysEye);
		SubSysCdag1.kronO(oldL.SubSysCdag1, oldR.SubSysEye);
	}

	//the rest is not needed in the superblock.

}




void Sub::trunc(const OP& truncU)
{

	SubSys.trunc(truncU);

	SubSysEye.trunc(truncU);

	SubSysC.trunc(truncU);

	SubSysCdag.trunc(truncU);





	if (SubSysC1.QDim.size() != 0)
	{
		SubSysC1.trunc(truncU);
		SubSysCdag1.trunc(truncU);
	}



}




void Sub::operator=(const Sub& block)
{
	clear();

	Orbital = block.Orbital;
	QorRl = block.QorRl;
	QorRr = block.QorRr;
	SubSys = block.SubSys;
	SubSysC = block.SubSysC;
	SubSysCdag = block.SubSysCdag;
	SubSysEye = block.SubSysEye;




	SubSysCdag1 = block.SubSysCdag1;
	SubSysC1 = block.SubSysC;


}

void Sub::show() const
{
	std::cout << "Orbital = " << Orbital << std::endl;
	std::cout << "QorRl = " << QorRl << std::endl;
	std::cout << "QorRr = " << QorRr << std::endl;
	std::cout << "SubSys : " << std::endl;
	SubSys.show();
	std::cout << "SubSysCdag : " << std::endl;
	SubSysCdag.show();
	std::cout << "SubSysC : " << std::endl;
	SubSysC.show();
	std::cout << "SubSysEye : " << std::endl;
	SubSysEye.show();








	std::cout << "SubSysCdag1: " << std::endl;
	SubSysCdag1.show();
	std::cout << "SubSysC1: " << std::endl;
	SubSysC1.show();




}





void Sub::clear()
{
	SubSys.clear();
	SubSysCdag.clear();
	SubSysC.clear();
	SubSysEye.clear();
	SubSysC1.clear();
	SubSysCdag1.clear();
}


//save the Sub.
void Sub::save()
{
	std::string str = itos(Orbital);
	str = "./data/" + str;
	std::ofstream outfile(str);
	outfile << Orbital << std::endl;
	outfile << QorRl << std::endl;
	outfile << QorRr << std::endl;
	SubSys.save(outfile);
	SubSysCdag.save(outfile);
	SubSysC.save(outfile);
	SubSysEye.save(outfile);
	SubSysC1.save(outfile);
	SubSysCdag1.save(outfile);
	outfile.close();

}




//read the Sub.
void Sub::read(int orbital_)
{
	clear();
	std::string str = itos(orbital_);
	str = "./data/" + str;
	std::ifstream infile(str);
	if (!infile.is_open())
	{
		std::cout << "the file " << str << " could not open. " << std::endl;
		exit(1);
	}
	infile >> Orbital;
	infile >> QorRl;
	infile >> QorRr;
	SubSys.read(infile);
	SubSysCdag.read(infile);
	SubSysC.read(infile);
	SubSysEye.read(infile);
	SubSysC1.read(infile);
	SubSysCdag1.read(infile);
	infile.close();
}




