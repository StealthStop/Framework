#ifndef MakeMT2Hemispheres_h
#define MakeMT2Hemispheres_h

#include "Framework/Framework/include/MT2Hemispheres.h"
#include "Framework/Framework/include/MT2HemispheresUtilities.h"
#include "Framework/Framework/include/Davismt2.h"
#include "Framework/Framework/include/TMctLib.h"
#include "Framework/Framework/include/mctlib.h"

#include <vector>
#include <iostream>
#include <cmath>


class MakeMT2Hemispheres
{

private:

    std::string myVarSuffix_;

    void getHemispheres(NTupleReader& tr) 
    {
        
    }

public:
    
    MakeMT2Hemispheres(std::string myVarSuffix = "")
        : myVarSuffix_       (myVarSuffix)
    {
        std::cout<<"Setting up MT2Hemispheres"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        getHemispheres(tr);
    }

};





#endif
