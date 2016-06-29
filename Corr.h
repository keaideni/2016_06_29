#include "Super.h"


#ifndef CORR_H
#define CORR_H


class Corr
{

private:
        
        

        int Orbital;
        int Type;//this is for the three types of operators.1 is for O, 2 is for ODag, 3 is for N.

public:
        OP CorrO;
        
        Corr();
        Corr(const Corr& corr_);

        void Initial(const Sub& sub, const Sub& subm, const int& orbital, const int& type, const OP& truncU);


        OP corro();
        

        int orbital();
        int type();


        void update(const Sub& m, const OP& truncU, const int& way);//way is for the site m on the block left or right.1 is for right and the other is for left.


        void read(const int& orbital, const int& type);
        void save();
        void show();
        void clear();





};














#endif