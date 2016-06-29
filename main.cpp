#include "Parameter.h"
#include "OP.h"
#include "SuperEnergy.h"
#include "DMRGP.h"



int OP::Max;
int main()
{
	Parameter para;
	para.read();
	
	OP::Max = para.ParticleNo + 1;

	
	//para.D = 200;
	





	DMRGP DMRG(para);


	
	system("pause");
}