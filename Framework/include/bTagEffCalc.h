//Tools Headers

#include "SusyAnaTools/Tools/SATException.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

//ROOT headers
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"
#include "TChain.h"

//STL Headers
#include<iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <sstream>
using namespace std;

const int nPtBins = 17;
const int nEtaBins = 3;

const double ptBins[]    =  {20,30,40,50,60,70,80,100,120,160,210,260,320,400,500,600,800,99999};
const double etaBins[]   =  {0.0,0.8,1.6,2.4};

const double csv_btag    = 0.8484; //Define here directly

const std::string spec = "bTagEff";

// === Main Function ===================================================
int main(int argc, char* argv[]) 
{
    TChain *fChain = 0;
    
    if (argc < 5) {
        std::cerr <<"Please give 4 arguments "<<"SubsampleName"<< " MaxEvent"<< "Startfile"<<" No. of Files to run" <<std::endl;
        std::cerr <<" Valid configurations are " << std::endl;
        std::cerr <<" ./remnantSys TTbarInc 1000 1 0" << std::endl;
        return -1;
	}

    const char *subSampleName = argv[1];
    const char *maxevent = argv[2];
    
    int numFiles = -1;
    int startFile = 0;
    // Change string arg to int
    int  maxEvent =  std::atoi(maxevent);
    numFiles =  std::atoi(argv[3]);
    startFile =  std::atoi(argv[4]);
    
    // Prepare file list and finalize it
    //TChain *fChain = 0;
    TString subSampleNameT = subSampleName;

    const std::string stFile = std::to_string(startFile);
    TString stFileT = stFile;
    TString specT = spec;

    TString FileName = subSampleNameT+"_"+spec+"_"+stFileT;
    std::cout<<FileName<<std::endl;
    TFile *f = nullptr;
    f =  new TFile(FileName+".root", "RECREATE");
      
    /*************************************************************/      
    //Declare efficiency histograms. n-> numerator d-> denominator
    // b-> bquark, c-> charm, usdg-> other light quarks.
    /*************************************************************/

    TH1::AddDirectory(kFALSE);
    TH2* n_eff_b =    new TH2D("n_eff_b", "bTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);
    TH2* n_eff_c =    new TH2D("n_eff_c", "cTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);
    TH2* n_eff_udsg = new TH2D("n_eff_udsg", "udsgTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);
    TH2* d_eff_b =    new TH2D("d_eff_b", "bTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);
    TH2* d_eff_c =    new TH2D("d_eff_c", "cTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);
    TH2* d_eff_udsg = new TH2D("d_eff_udsg", "udsgTag_Efficiency"+stFileT, nPtBins, ptBins, nEtaBins, etaBins);

    // TH2* n_eff_b =    new TH2D("n_eff_b_"+samplesT, "bTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
    // TH2* n_eff_c =    new TH2D("n_eff_c_"+samplesT, "cTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
    // TH2* n_eff_udsg = new TH2D("n_eff_udsg_"+samplesT, "udsgTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
    // TH2* d_eff_b =    new TH2D("d_eff_b_"+samplesT, "bTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
    // TH2* d_eff_c =    new TH2D("d_eff_c_"+samplesT, "cTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
    // TH2* d_eff_udsg = new TH2D("d_eff_udsg_"+samplesT, "udsgTag_Efficiency"+samplesT, nPtBins, ptBins, nEtaBins, etaBins);
      
    n_eff_b->GetXaxis()->SetTitle("p_{T} [GeV]");
    n_eff_b->GetYaxis()->SetTitle("#eta");

    n_eff_c->GetXaxis()->SetTitle("p_{T} [GeV]");
    n_eff_c->GetYaxis()->SetTitle("#eta");
      
    n_eff_udsg->GetXaxis()->SetTitle("p_{T} [GeV]");
    n_eff_udsg->GetYaxis()->SetTitle("#eta");
 

    //const string condor =  (argc == 6) ? argv[5]: "";

    //condor.empty() ? AnaSamples::SampleSet():AnaSamples::SampleSet(argv[5]);  

    //AnaSamples::SampleSet        ss("sampleSets.txt");
    //AnaSamples::SampleCollection sc("sampleCollections.txt", ss);
        /*                           
    double ScaleMC = 1.;                                                                              
    if(ss[subSampleName] != ss.null()) {
        
        fChain = new TChain(ss[subSampleName].treePath.c_str());
        ss[subSampleName].addFilesToChain(fChain, startFile, numFiles);

	    ScaleMC = ss[subSampleName].getWeight();
    } */

    //AnaFunctions::prepareForNtupleReader();
    //NTupleReader *tr =0;
    //tr = new NTupleReader(fChain);
    
    // --- Analyse events --------------------------------------------
    std::cout<<"First loop begin: "<<std::endl;
    int entries = tr->getNEntries();
    std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<ScaleMC<<std::endl; 
    std::cout<<"maxevent: "<<maxEvent<<std::endl;

      
    /*************************************************************/
    // Event loop begins                                                                                                            
    /*************************************************************/
/*    while(tr->getNextEvent())
    {

	if(maxEvent>=0 && tr->getEvtNum() > maxEvent ) break;
	// Add print out of the progress of looping
	if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) 
	  {
	    std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
	  }
	  
	  
        const vector<TLorentzVector> inputJets = tr->getVec<TLorentzVector>("Jets");
        const vector<double> recoJetsBtag = tr->getVec<double>("Jets_bDiscriminatorCSV");
        const vector<int> recoJetsFlavor = tr->getVec<int>("Jets_hadronFlavor");
         
        double iniWeight = tr->getVar<double>("evtWeight");

        double stored_weight = subSampleNameT.Contains("Data") ? 1 : tr->getVar<double>("stored_weight");
        int sign_of_stored_weight = (stored_weight > 0) ? 1 : ((stored_weight < 0) ? -1 : 0);

        double evtWeight = iniWeight >=0 ? iniWeight * sign_of_stored_weight : iniWeight;
     
        for(unsigned int ij=0; ij<inputJets.size(); ij++)
        {
            double pt = inputJets[ij].Pt();
            double eta = fabs(inputJets[ij].Eta());
            if( ! AnaFunctions::jetPassCuts(inputJets[ij], AnaConsts::bTagArr) ) continue;
            double csv = recoJetsBtag.at(ij);
            int flav =  abs(recoJetsFlavor.at(ij));
	      
            if(flav==5) //b Jets
            {
                d_eff_b->Fill(pt,eta, evtWeight);
                if(csv > csv_btag) n_eff_b->Fill(pt,eta, evtWeight);
            }
            else if(flav==4) // c jets
            {
                d_eff_c->Fill(pt,eta, evtWeight);
                if(csv > csv_btag) n_eff_c->Fill(pt,eta, evtWeight);
            }
            else if(flav<4 || flav==21) // other flavours
            {
                d_eff_udsg->Fill(pt,eta, evtWeight);
                if(csv > csv_btag) n_eff_udsg->Fill(pt,eta, evtWeight);
            }
        }
    }
*/
    /*************************************************************/
    // End of Event loop
    /*************************************************************/
/*  
    f->cd();
    d_eff_b->Write();
    d_eff_c->Write();
    d_eff_udsg->Write();


    n_eff_b->Write();
    n_eff_c->Write();
    n_eff_udsg->Write();


    f->Close();


  fChain->Reset();
  return 0;
*/
 
}

/**************************************************************************/
// End of Main function
// MakeFile example to compile this is at 
// https://github.com/humkies/bTag/blob/master/Makefile#L67
/**************************************************************************/





 





