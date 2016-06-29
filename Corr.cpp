#include "Corr.h"

Corr::Corr(){}

Corr::Corr(const Corr& corr)
{
        Orbital = corr.Orbital;

        CorrO = corr.CorrO;

        Type = corr.Type;

}




void Corr::Initial(const Sub& sub,const Sub& subm, const int& orbital, const int& type, const OP& truncU_)
{
        clear();
        Orbital = orbital;
        Type = type;
        OP tempCorrO;
        if(Type == 1)
        {
                tempCorrO = subm.SubSysC;

        }else if(Type == 2)
        {
                tempCorrO = subm.SubSysCdag;
        }else
        {
                tempCorrO.time(subm.SubSysCdag, subm.SubSysC);
        }

        CorrO.kronO(sub.SubSysEye, tempCorrO);
        CorrO.trunc(truncU_);


}




OP Corr::corro()
{
        return CorrO;
}




int Corr::orbital()
{
        return Orbital;
}



int Corr::type()
{
        return Type;
}



void Corr::update(const Sub& m, const OP& truncU, const int& way)
{
        if(way == 1)
        {
               OP tempCorrO;
               tempCorrO.kronO(this->CorrO, m.SubSysEye);
               tempCorrO.trunc(truncU); 
               CorrO = tempCorrO;
        }else
        {
                OP tempCorrO;
                tempCorrO.kronO(m.SubSysEye, this->CorrO);

                tempCorrO.trunc(truncU); 
                CorrO = tempCorrO;
        }
}





void Corr::save()
{
        std::string str = itos(Orbital);
        if(Type == 1)
        {
                str = "./Corr/C" + str;
        }else if(Type == 2)
        {
                str = "./Corr/CDag" + str;
        }else
        {
                str = "./Corr/N" + str;
        }

        std::ofstream outfile(str);
        
        outfile << Orbital <<std::endl;
        outfile << Type <<std::endl;
        CorrO.save(outfile);

        outfile.close();
}



void Corr::read(const int& orbital, const int& type)
{
        clear();
        std::string str = itos(orbital);
        if(type == 1)
        {
                str = "./Corr/C" + str;
        }else if(type == 2)
        {
                str = "./Corr/CDag" + str;
        }else
        {
                str = "./Corr/N" + str;
        }

        std::ifstream infile(str);
        if(!infile.is_open())
        {
                std::cout<<"the file "<<str<<" doesn't exist!"<<std::endl;
                exit(1);
        }
        infile >> Orbital;
        infile >> Type;

        CorrO.read(infile);


}




void Corr::show()
{
        std::string str = itos(Orbital);
        if(Type == 1)
        {
                str = "the Operator C on the msite "+str;
        }else if(Type == 2)
        {
                str = "the Operator CDag on the msite "+str;
        }else
        {
                str = "the Operator N on the msite "+str;
        }
        std::cout << str <<std::endl;
        CorrO.show();

}



void Corr::clear()
{
        CorrO.clear();
}




