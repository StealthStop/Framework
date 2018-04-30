#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TList.h"
#include "TLegendEntry.h"
#include "string.h"
#include <iostream>
#include "TTreeReader.h"
#include "TChain.h"
#include "TProfile.h"
#include <cmath>
#include "TProfile2D.h"
#include "TH2D.h"

void CombineTree(const std::string& infile, const std::string& outfile)
{ 
    gSystem->Exec( ("rm  " + outfile).c_str() );
   
    std::ifstream in;
    in.open( infile.c_str() );
    if(!in.is_open())
    {
        std::cout << "Cannot open list file: " << infile << std::endl;
        return;  
    }
    
    TChain* chain = new TChain("mvatraintt");
    
    std::string line;
    while( in.good() )
    {
        if( !std::getline(in,line) ) break; // We read a line from the file
        if( !chain->Add(line.c_str()) )
        {
            std::cout << "Problem loading tree from " << line << std::endl;
        }
        else
        {
            std::cout << "Adding file: " << line << "..." << std::endl;
        }
        
    }
    
    in.close();
    
    chain->Merge( outfile.c_str() );
    
    delete chain;
}
