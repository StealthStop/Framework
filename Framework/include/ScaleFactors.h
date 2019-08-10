#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
#include "TFile.h"
#include "TKey.h"
#include "Framework/Framework/include/Utility.h"

class ScaleFactors
{
private:
    std::string myVarSuffix_;

    std::shared_ptr<TH2F> eleSFHistoTight_;        
    std::shared_ptr<TH2F> eleSFHistoIso_;
    std::shared_ptr<TH2F> eleSFHistoIP2D_;
    std::shared_ptr<TH2F> eleSFHistoReco_;
    std::shared_ptr<TH2F> eleSFHistoTrig_;
    std::shared_ptr<TH2F> muSFHistoMedium_;
    std::shared_ptr<TH2F> muSFHistoIso_;
    std::shared_ptr<TH2F> muSFHistoTrig_;
    std::shared_ptr<TGraph> muSFHistoReco_;
    std::shared_ptr<TH1F> puSFHisto_;
    std::shared_ptr<TH1F> puSFUpHisto_;
    std::shared_ptr<TH1F> puSFDownHisto_;
    std::shared_ptr<TH2F> L1Prefireing_;
    std::map<std::string, double> sfMeanMap_;

    template<typename Th, typename Tb> const int findBin(const std::shared_ptr<Th>& h, const Tb v, const std::string& axis)
    {
        int bin = -1;
        if(axis=="X")
        {
            bin = h->GetXaxis()->FindBin(v);
            if( v >= h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()) ) bin = h->GetNbinsX();
        }
        else if(axis=="Y")
        {
            bin = h->GetYaxis()->FindBin(v);
            if( v >= h->GetYaxis()->GetBinUpEdge(h->GetNbinsY()) ) bin = h->GetNbinsY();           
        }
        if(bin == -1) std::cerr<<"There was an error in extracting the bin index for a scale factor"<<std::endl;
        return bin;
    }

    void scaleFactors(NTupleReader& tr)
    {
        // --------------------------------------------------------------------------------------
        // Calculate scale SF and variation
        // Following the example in SusyAnaTools PDFUncertainty.h
        // The scale weights are calculated using the envelope method and we ignore all anti-correlated variations (5 and 7)
        // --------------------------------------------------------------------------------------
        const auto& scaleWeights = tr.getVec<double>("ScaleWeights");
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& runYear      = tr.getVar<std::string>("runYear");
        
        double scaleWeightNominal = 1.0;
        std::vector<double>  myScaleWeights(6, 1.0);
        if( scaleWeights.size() == 9 ) 
        {
            //If there are not exactly 9 scale factors, then the vector  was filled incorrectly
            scaleWeightNominal = scaleWeights[0];
            myScaleWeights[0]  = scaleWeights[1];
            myScaleWeights[1]  = scaleWeights[2];
            myScaleWeights[2]  = scaleWeights[3];
            myScaleWeights[3]  = scaleWeights[4];
            myScaleWeights[4]  = scaleWeights[6];
            myScaleWeights[5]  = scaleWeights[8];
        }
        
        double scaleWeightUpperBound = *std::max_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        double scaleWeightLowerBound = *std::min_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        if( !std::isfinite(scaleWeightNominal) || !std::isfinite(scaleWeightUpperBound) || !std::isfinite(scaleWeightLowerBound) ) 
        {
            //TODO: There are some NaN/Inf values in the Diboson channel - still need to figure out why this is an issue.
            scaleWeightNominal = 1.0;
            scaleWeightUpperBound = 1.0;
            scaleWeightLowerBound = 1.0;
        }

        double scaleWeightUpperBound_corr = scaleWeightUpperBound;
        double scaleWeightLowerBound_corr = scaleWeightLowerBound;
        if(sfMeanMap_.find(filetag+"_sclUp") != sfMeanMap_.end() && sfMeanMap_.find(filetag+"_sclDown") != sfMeanMap_.end()) 
        {            
            const double meanUp = sfMeanMap_[filetag+"_sclUp"];
            const double meanDown = sfMeanMap_[filetag+"_sclDown"];
            scaleWeightUpperBound_corr = (1/meanUp)*scaleWeightUpperBound;
            scaleWeightLowerBound_corr = (1/meanDown)*scaleWeightLowerBound;
        }

        tr.registerDerivedVar("scaleWeightUp"+myVarSuffix_,   scaleWeightUpperBound_corr);
        tr.registerDerivedVar("scaleWeightDown"+myVarSuffix_, scaleWeightLowerBound_corr);
        tr.registerDerivedVar("scaleWeightUpUncor"+myVarSuffix_,   scaleWeightUpperBound);
        tr.registerDerivedVar("scaleWeightDownUncor"+myVarSuffix_, scaleWeightLowerBound);
        tr.registerDerivedVar("scaleWeightNom"+myVarSuffix_,  scaleWeightNominal);

        // --------------------------------------------------------------------------------------
        // Calculate parton shower variation
        // Note: not all samples have these weights stored, give them default value of 1. 
        // --------------------------------------------------------------------------------------        
        double PSweight_ISRUp = 1.0, PSweight_ISRDown = 1.0; 
        double PSweight_FSRUp = 1.0, PSweight_FSRDown = 1.0; 
        double PSweight_ISRUp_2 = 1.0, PSweight_ISRDown_2 = 1.0; 
        double PSweight_FSRUp_2 = 1.0, PSweight_FSRDown_2 = 1.0; 
        if(tr.hasVar("PSweights"))
        {
            const auto& PSweights = tr.getVec<double>("PSweights");
            if(PSweights.size() >= 12) // should have size of 14, but just put 12 or more to be able to use the sample with the bug
            {
                // Get nominal one so we can normalize it
                double MEweight = PSweights[0];
                // reduced variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/sqrt(2) and sqrt(2)
                PSweight_ISRUp = PSweights[2]/MEweight;
                PSweight_FSRUp = (PSweights[3]/MEweight < 10.0) ? PSweights[3]/MEweight : 1.0;
                PSweight_ISRDown = PSweights[4]/MEweight;
                PSweight_FSRDown = (PSweights[5]/MEweight < 10.0) ? PSweights[5]/MEweight : 1.0;
                // nominal variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/2 and 2
                PSweight_ISRUp_2 = PSweights[6]/MEweight;
                PSweight_FSRUp_2 = (PSweights[7]/MEweight < 10.0) ? PSweights[7]/MEweight : 1.0;
                PSweight_ISRDown_2 = PSweights[8]/MEweight;
                PSweight_FSRDown_2 = (PSweights[9]/MEweight < 10.0) ? PSweights[9]/MEweight : 1.0;
            }
        }
        tr.registerDerivedVar("PSweight_ISRUp"+myVarSuffix_,   PSweight_ISRUp);
        tr.registerDerivedVar("PSweight_ISRDown"+myVarSuffix_, PSweight_ISRDown);
        tr.registerDerivedVar("PSweight_FSRUp"+myVarSuffix_,   PSweight_FSRUp);
        tr.registerDerivedVar("PSweight_FSRDown"+myVarSuffix_, PSweight_FSRDown);
        tr.registerDerivedVar("PSweight_ISRUp_2"+myVarSuffix_,   PSweight_ISRUp_2);
        tr.registerDerivedVar("PSweight_ISRDown_2"+myVarSuffix_, PSweight_ISRDown_2);
        tr.registerDerivedVar("PSweight_FSRUp_2"+myVarSuffix_,   PSweight_FSRUp_2);
        tr.registerDerivedVar("PSweight_FSRDown_2"+myVarSuffix_, PSweight_FSRDown_2);

        // --------------------------------------------------------------------------------------
        // Now calculate the PDF scale factor and uncertainty 
        // Based on the 100 different replica values stored in PDFweights using envelope method and the median
        // --------------------------------------------------------------------------------------
        const auto& PDFweights = tr.getVec<double>("PDFweights");
        double central = 1.0, NNPDF_from_median_up = 1.0, NNPDF_from_median_down = 1.0;
        if(PDFweights.size() > 0)
        {
            const double reqCL = 0.68; //Choose a confidence level for the uncertainty
            std::vector<double> sortedPDFWeights = PDFweights; //Cannot sort a constant
            std::sort( sortedPDFWeights.begin() + 1, sortedPDFWeights.end() );
        
            const int upper = std::round( 0.5*(1 + reqCL)*100.0 );
            const int lower = 1 + std::round( 0.5*(1 - reqCL)*100.0 );

            central  = 0.5*( sortedPDFWeights[50] + sortedPDFWeights[51] ); //Exactly 100 entries
            const double errminus = abs(central - sortedPDFWeights[lower]);
            const double errplus  = abs(central - sortedPDFWeights[upper]);
            double errsymm  = 0.5*( errplus + errminus );

            NNPDF_from_median_up   = central + errplus;
            NNPDF_from_median_down = central - errplus;
            NNPDF_from_median_up   = NNPDF_from_median_up/central   > 2.0 ? 1.0 : NNPDF_from_median_up/central   < -2.0 ? 1.0 : NNPDF_from_median_up/central;        
            NNPDF_from_median_down = NNPDF_from_median_down/central > 2.0 ? 1.0 : NNPDF_from_median_down/central < -2.0 ? 1.0 : NNPDF_from_median_down/central;
        }
        if( !std::isfinite(central) || !std::isfinite(NNPDF_from_median_up) || !std::isfinite(NNPDF_from_median_down) ) 
        {
            NNPDF_from_median_up    = 1.0;
            NNPDF_from_median_down  = 1.0;
            central                 = 1.0;
        }
        
        double NNPDF_from_median_up_corr = NNPDF_from_median_up;
        double NNPDF_from_median_down_corr = NNPDF_from_median_down;
        if(sfMeanMap_.find(filetag+"_pdf_Up") != sfMeanMap_.end() && sfMeanMap_.find(filetag+"_pdf_Down") != sfMeanMap_.end()) 
        {
            const double meanUp = sfMeanMap_[filetag+"_pdf_Up"];
            const double meanDown = sfMeanMap_[filetag+"_pdf_Down"];
            NNPDF_from_median_up_corr = (1/meanUp)*NNPDF_from_median_up;
            NNPDF_from_median_down_corr = (1/meanDown)*NNPDF_from_median_down;
        }

        tr.registerDerivedVar( "PDFweightUp"+myVarSuffix_,   NNPDF_from_median_up_corr );
        tr.registerDerivedVar( "PDFweightDown"+myVarSuffix_, NNPDF_from_median_down_corr );
        tr.registerDerivedVar( "PDFweightUpUncor"+myVarSuffix_,   NNPDF_from_median_up );
        tr.registerDerivedVar( "PDFweightDownUncor"+myVarSuffix_, NNPDF_from_median_down );
        tr.registerDerivedVar( "PDFweightNom"+myVarSuffix_,  central );
        
        // --------------------------------------------------------------------------------------
        // Now calculate the electron scale factor and uncertainty 
        // --------------------------------------------------------------------------------------
        const auto& electrons         = tr.getVec<TLorentzVector>("Electrons");
        const auto& goodElectrons     = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        
        double totGoodElectronSF      = 1.0, totGoodElectronSF_Up   = 1.0, totGoodElectronSF_Down = 1.0;
        double totGoodElectronSFErr   = 0.0, totGoodElectronSFPErr2 = 0.0;
        double noTrigGoodElectronSF   = 1.0, noTrigGoodElectronSFErr = 0.0, noTrigGoodElectronSFPErr2 = 0.0;
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) 
        {
            if( !goodElectrons.at(iel) ) continue;
            
            //Get the scale factor from the rootfile
            const double elpt  = electrons.at(iel).Pt();
            const double eleta = electrons.at(iel).Eta();
            int xbinElTight = -1, ybinElTight = -1, xbinElIso = -1, ybinElIso = -1, xbinElReco = -1, ybinElReco = -1, xbinElTrig = -1, ybinElTrig = -1;

            if( runYear == "2016") 
            {
                //Find the bin indices (binned by x: eta and y: pt ) for the 2016 scale factor for Tight ID ( has 30 bins total in 2D parameter space )
                //This is the same index as for the Iso ID
                xbinElTight = findBin<TH2F, double>(eleSFHistoTight_, eleta, "X");
                ybinElTight = findBin<TH2F, double>(eleSFHistoTight_, elpt,  "Y");
                //Since the binning for MiniIso < 0.1 is the same as that for Tight ID, we will use the same values (values initialized for uniformity later on). 
                //Same for the IP2D (hence you do not need an extra set of variables).
                xbinElIso = xbinElTight;
                ybinElIso = ybinElTight;
                //Find the bin indices (binned by x: eta and y: pt ) for the 2016 scale factor for Data/MC differences (reco eff) ( has 90 bins total in 2D parameter space )
                xbinElReco = findBin<TH2F, double>(eleSFHistoReco_, eleta, "X");
                ybinElReco = findBin<TH2F, double>(eleSFHistoReco_, elpt,  "Y");
                //Find the bin indices (binned by x: pt and y: eta ) for the 2016 scale factor for Data/MC trig efficiency
                xbinElTrig = findBin<TH2F, double>(eleSFHistoTrig_, elpt,  "X");
                ybinElTrig = findBin<TH2F, double>(eleSFHistoTrig_, eleta, "Y");
            }
            else if( runYear == "2017" ) 
            {                
                //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for Tight ID ( has 180 bins total in 2D parameter space )
                xbinElTight = findBin<TH2F, double>(eleSFHistoTight_, eleta, "X");
                ybinElTight = findBin<TH2F, double>(eleSFHistoTight_, elpt,  "Y");
                //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for MiniIso of 0.1  ( has 210 bins total in 2D parameter space )
                xbinElIso = findBin<TH2F, double>(eleSFHistoIso_, eleta, "X");
                ybinElIso = findBin<TH2F, double>(eleSFHistoIso_, elpt,  "Y");
                //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for Data/MC comparison (reco eff) ( has 144 bins total in 2D parameter space )
                xbinElReco = findBin<TH2F, double>(eleSFHistoReco_, eleta, "X");
                ybinElReco = findBin<TH2F, double>(eleSFHistoReco_, elpt,  "Y");
                //Find the bin indices (binned by x: pt and y: eta ) for the 2017 scale factor for Data/MC trig efficiency
                xbinElTrig = findBin<TH2F, double>(eleSFHistoTrig_, elpt,  "X");
                ybinElTrig = findBin<TH2F, double>(eleSFHistoTrig_, eleta, "Y");
            }
            else if( runYear == "2018" )
            {
                //Find the bin indices (binned by x: eta and y: pt ) for the 2018 scale factor for Tight ID ( has 96 bins total in 2D parameter space )
                xbinElTight = findBin<TH2F, double>(eleSFHistoTight_, eleta, "X");
                ybinElTight = findBin<TH2F, double>(eleSFHistoTight_, elpt,  "Y");                
                //Find the bin indices (binned by x: eta and y: pt ) for the 2018 scale factor for MiniIso of 0.1  ( has 96 bins total in 2D parameter space )
                xbinElIso = findBin<TH2F, double>(eleSFHistoIso_, eleta, "X");
                ybinElIso = findBin<TH2F, double>(eleSFHistoIso_, elpt,  "Y");
                //Find the bin indices (binned by x: eta and y: pt ) for the 2018 scale factor for Data/MC comparison (reco eff) ( has 98 bins total in 2D parameter space )
                xbinElReco = findBin<TH2F, double>(eleSFHistoReco_, eleta, "X");
                ybinElReco = findBin<TH2F, double>(eleSFHistoReco_, elpt,  "Y");
            }

            if( xbinElTight != -1 && ybinElTight != -1 && xbinElIso != -1 && ybinElIso != -1 && xbinElReco != -1 && ybinElReco != -1 ) 
            {
                const double eleTightSF       = eleSFHistoTight_->GetBinContent( xbinElTight, ybinElTight );
                const double eleTightSFErr    = eleSFHistoTight_->GetBinError( xbinElTight, ybinElTight );
                const double eleTightPErr     = eleTightSFErr/eleTightSF;

                double eleIsoSF         = eleSFHistoIso_->GetBinContent( xbinElIso, ybinElIso );
                double eleIsoSFErr      = eleSFHistoIso_->GetBinError( xbinElIso, ybinElIso );
                double eleIsoPErr       = eleIsoSFErr/eleIsoSF;

                const double eleRecoSF        = eleSFHistoReco_->GetBinContent( xbinElReco, ybinElReco );
                const double eleRecoSFErr     = eleSFHistoReco_->GetBinError( xbinElReco, ybinElReco );
                const double eleRecoPErr      = eleRecoSFErr/eleRecoSF;

                const double eleTrigSF        = eleSFHistoTrig_->GetBinContent( xbinElTrig, ybinElTrig );
                const double eleTrigSFErr     = eleSFHistoTrig_->GetBinError( xbinElTrig, ybinElTrig );
                const double eleTrigPErr      = eleTrigSFErr/eleTrigSF;
                
                if( runYear == "2016" ) 
                { 
                    //The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    //If this is the year 2016, we need to add the IP2D histogram scale factors into the Iso scale factor
                    const double eleIP2DSF    = eleSFHistoIP2D_->GetBinContent( xbinElIso, ybinElIso );
                    const double eleIP2DSFErr = eleSFHistoIP2D_->GetBinError( xbinElIso, ybinElIso );
                    const double eleIP2DPErr  = eleIP2DSFErr/eleIP2DSF;

                    eleIsoSF            = eleIsoSF*eleIP2DSF;
                    eleIsoPErr          = utility::addInQuad( eleIsoPErr, eleIP2DPErr );
                    eleIsoSFErr         = eleIsoPErr*eleIsoSF;
                }
                
                const double eleNoTrigSF      = eleTightSF*eleIsoSF*eleRecoSF;
                const double eleTotSF         = eleNoTrigSF*eleTrigSF; 

                const double eleNoTrigPErr    = utility::addInQuad( eleTightPErr, eleIsoPErr, eleRecoPErr );
                const double eleTotPErr       = utility::addInQuad( eleNoTrigPErr, eleTrigPErr );

                const double eleTotSFErr      = eleTotPErr*eleTotSF;
                const double eleNoTrigSFErr   = eleNoTrigPErr*eleNoTrigSF;

                totGoodElectronSF       *= eleTotSF; 
                noTrigGoodElectronSF    *= eleNoTrigSF;

                totGoodElectronSFPErr2  += eleTotPErr*eleTotPErr;
                noTrigGoodElectronSFPErr2 += eleNoTrigPErr*eleNoTrigPErr;
            }            
        }

        totGoodElectronSFErr    = std::sqrt(totGoodElectronSFPErr2) * totGoodElectronSF;
        noTrigGoodElectronSFErr = std::sqrt(noTrigGoodElectronSFPErr2) * noTrigGoodElectronSF;

        totGoodElectronSF_Up    = totGoodElectronSF + totGoodElectronSFErr;
        totGoodElectronSF_Down  = totGoodElectronSF - totGoodElectronSFErr;

        tr.registerDerivedVar( "totGoodElectronSF"+myVarSuffix_,      totGoodElectronSF );
        tr.registerDerivedVar( "totGoodElectronSFErr"+myVarSuffix_,   totGoodElectronSFErr );
        tr.registerDerivedVar( "totGoodElectronSF_Up"+myVarSuffix_,   totGoodElectronSF_Up );
        tr.registerDerivedVar( "totGoodElectronSF_Down"+myVarSuffix_, totGoodElectronSF_Down );

        tr.registerDerivedVar( "noTrigGoodElectronSF"+myVarSuffix_,   noTrigGoodElectronSF );
        tr.registerDerivedVar( "noTrigGoodElectronSFErr"+myVarSuffix_, noTrigGoodElectronSFErr );

        // --------------------------------------------------------------------------------------
        // Adding code for implementing muon scale factors
        // --------------------------------------------------------------------------------------
        const auto& muons         = tr.getVec<TLorentzVector>("Muons");
        const auto& goodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_);

        double totGoodMuonSF      = 1.0, totGoodMuonSF_Up   = 1.0, totGoodMuonSF_Down = 1.0;
        double totGoodMuonSFErr   = 0.0, totGoodMuonSFPErr2 = 0.0;

        double noTrigGoodMuonSF   = 1.0, noTrigGoodMuonSFErr = 0.0, noTrigGoodMuonSFPErr2 = 0.0;
        for( unsigned int imu = 0; imu < muons.size(); imu++ ) 
        {            
            if( !goodMuons.at(imu) ) continue;
            
            //Get the scale factor from the rootfile
            const double mupt = muons.at(imu).Pt();
            const double mueta = muons.at(imu).Eta();
            int xbinMuMedium = -1, ybinMuMedium = -1, xbinMuIso = -1, ybinMuIso = -1, xbinMuTrig = -1, ybinMuTrig = -1;

            if( runYear == "2016" ) 
            {                
                //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2016 scale factor for both Medium ID and MiniIso of 0.2 (has only 20 bins) - same binning!
                xbinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, mupt,  "X");
                ybinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, abs(mueta), "Y");
                //For consistency with 2017, copy the values over to these variables.
                xbinMuIso = xbinMuMedium;
                ybinMuIso = ybinMuMedium;
                //Find the bin indices (binned by x: pt and y: eta ) for the 2016 scale factor for Data/MC Trigger Efficiency
                xbinMuTrig = findBin<TH2F, double>(muSFHistoTrig_, mupt,  "X");
                ybinMuTrig = findBin<TH2F, double>(muSFHistoTrig_, mueta, "Y");
            } 
            else if( runYear == "2017" ) 
            {
                //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2017 scale factor for Medium ID (has only 24 bins)
                xbinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, mupt,       "X");
                ybinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, abs(mueta), "Y");
                //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2017 scale factor for MiniIso of 0.2 (has only 20 bins)
                xbinMuIso = findBin<TH2F, double>(muSFHistoIso_, mupt,       "X");
                ybinMuIso = findBin<TH2F, double>(muSFHistoIso_, abs(mueta), "Y");
                //Find the bin indices (binned by x: pt and y: eta ) for the 2017 scale factor for Data/MC Trigger Efficiency
                xbinMuTrig = findBin<TH2F, double>(muSFHistoTrig_, mupt,  "X");
                ybinMuTrig = findBin<TH2F, double>(muSFHistoTrig_, mueta, "Y");
            }
            else if( runYear == "2018" )
            {
                //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2018 scale factor for Medium ID (has only 36 bins)
                xbinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, mupt,       "X");
                ybinMuMedium = findBin<TH2F, double>(muSFHistoMedium_, abs(mueta), "Y");
                //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2018 scale factor for MiniIso of 0.2 (has only 20 bins)
                //The muon POG will not release MiniIso scale factors until the UltraLegacy, so we recommend to use the 2017 Data/FullSim SFs for MiniIso also for 2018
                xbinMuIso = findBin<TH2F, double>(muSFHistoIso_, mupt,       "X");
                ybinMuIso = findBin<TH2F, double>(muSFHistoIso_, abs(mueta), "Y");            
            }

            if( xbinMuMedium != 0 && ybinMuMedium != 0 && xbinMuIso != 0 && ybinMuIso != 0 ) 
            {
                //The SUSLepton Twiki claims that the errors in the histogrm are purely statistical and can be ignored and recommends a 3% error for each leg (ID+IP+ISO)
                const double muMediumSF     = muSFHistoMedium_->GetBinContent( xbinMuMedium, ybinMuMedium );
                const double muIsoSF        = muSFHistoIso_->GetBinContent( xbinMuIso, ybinMuIso );
                const double muTrigSF       = muSFHistoTrig_->GetBinContent( xbinMuTrig, ybinMuTrig );
                const double muTrigSFErr    = muSFHistoTrig_->GetBinError( xbinMuTrig, ybinMuTrig );
                const double muTrigSFPErr   = muTrigSFErr/muTrigSF;

                double muNoTrigSF     = muMediumSF * muIsoSF; 
                double muTotSF        = muNoTrigSF * muTrigSF;

                const double muNoTrigSFPErr2 = 0.03*0.03;
                const double muTotSFPErr2   = muNoTrigSFPErr2 + muTrigSFPErr*muTrigSFPErr;
                
                if( runYear == "2016" ) 
                {
                    //For the general track reconstruction they claim that the errors for the systematic still need to be finalized - does not seem to have been finalized as of Dec 2018
                    //This reconstruction value only exists for 2016 - SUS SF people say the 3% will include the reco scale factor uncertainty for now
                    const double muRecoSF   = muSFHistoReco_->Eval( mueta );
                    muTotSF           = muTotSF*muRecoSF;
                    muNoTrigSF        = muNoTrigSF*muRecoSF;
                }
                totGoodMuonSF      *= muTotSF;
                noTrigGoodMuonSF   *= muNoTrigSF;

                totGoodMuonSFPErr2 += muTotSFPErr2;
                noTrigGoodMuonSFPErr2 += muNoTrigSFPErr2;
            }
        }

        totGoodMuonSFErr    = std::sqrt(totGoodMuonSFPErr2)*totGoodMuonSF;
        totGoodMuonSF_Up    = totGoodMuonSF + totGoodMuonSFErr;
        totGoodMuonSF_Down  = totGoodMuonSF - totGoodMuonSFErr;

        noTrigGoodMuonSFErr = std::sqrt(noTrigGoodMuonSFPErr2)*noTrigGoodMuonSF;

        tr.registerDerivedVar( "totGoodMuonSF"+myVarSuffix_,      totGoodMuonSF );
        tr.registerDerivedVar( "totGoodMuonSFErr"+myVarSuffix_,   totGoodMuonSFErr );
        tr.registerDerivedVar( "totGoodMuonSF_Up"+myVarSuffix_,   totGoodMuonSF_Up );
        tr.registerDerivedVar( "totGoodMuonSF_Down"+myVarSuffix_, totGoodMuonSF_Down );

        tr.registerDerivedVar( "noTrigGoodMuonSF"+myVarSuffix_,   noTrigGoodMuonSF );
        tr.registerDerivedVar( "noTrigGoodMuonSFErr"+myVarSuffix_,noTrigGoodMuonSFErr );

        // --------------------------------------------------------------------------------------
        // Adding a scale factor that corrects the disagreement between data and MC for Ht
        // --------------------------------------------------------------------------------------
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& isSignal = tr.getVar<bool>("isSignal");

        auto htScaleFactor = [](int nJets, double HT, const std::string& runYear) 
        { 
            double norm = 0.0;
            double expo = 0.0;
            if( runYear == "2016" ) 
            {
                norm = 0.03422*nJets + 0.9367;
                expo = (-0.02310*nJets - 0.0940 )/1000;
            }
            else if( runYear == "2017" ) 
            {
                norm = 0.02565*nJets + 0.9635;
                expo = (-0.01418*nJets - 0.1101 )/1000;
            }
            return norm*exp( expo*HT ); 
        };
        auto htScaleFactorFlat2000 = [](int nJets, double HT, const std::string& runYear)
        {
            double norm = 0.0;
            double expo = 0.0;
            if( runYear == "2016" ) 
            {
                norm = 0.03422*nJets + 0.9367;
                expo = (-0.02310*nJets - 0.0940 )/1000;
            }
            else if( runYear == "2017" ) 
            {
                norm = 0.02565*nJets + 0.9635;
                expo = (-0.01418*nJets - 0.1101 )/1000;
            }
            
            if( HT > 2000 )
                return norm*exp( expo*2000.00 );
            else
                return norm*exp( expo*HT ); 
        };
        auto htScaleFactorNJet7 = [](double HT, const std::string& runYear)
        {
            double norm = 0.0;
            double expo = 0.0;
            
            if( runYear == "2016" ) 
            {
                norm = 0.03422*7 + 0.9367;
                expo = (-0.02310*7 - 0.0940 )/1000;
            }
            else if( runYear == "2017" ) 
            {
                norm = 0.02565*7 + 0.9635;
                expo = (-0.01418*7 - 0.1101 )/1000;
            }            
            return norm*exp( expo*HT ); 
        };
        auto htScaleFactorMG = [](int nJets, double HT, std::string filetag)
        {
            double norm     = 0.0;
            double expo     = 0.0;

            if( filetag.find("2016") != std::string::npos ) {
                norm = 0.01802*nJets + 0.9762;
                expo = (-0.003885*nJets - 0.1074)/1000;
            }
            else {
                norm = 0.01818*nJets + 1.0535;
                expo = (0.00425*nJets - 0.3170)/1000;
            }

            return norm*exp( expo*HT );
        };

        double htDerivedweight = 1.0, htDerivedweightFlat2000 = 1.0, htDerivedweightNJet7 = 1.0, htDerivedweightMG    = 1.0;
        double htDerivedweightUncor = htScaleFactor(NGoodJets_pt30, HT_trigger_pt30, runYear);
        double htDerivedweightFlat2000Uncor = htScaleFactorFlat2000(NGoodJets_pt30, HT_trigger_pt30, runYear);
        double htDerivedweightNJet7Uncor = htScaleFactorNJet7(HT_trigger_pt30, runYear);
        double htDerivedweightMGUncor = htScaleFactorMG(NGoodJets_pt30, HT_trigger_pt30, filetag);

        double htScaleUp = 1.0, htScaleDown = 1.0, htScaleUpMG = 1.0, htScaleDownMG = 1.0;
        if( sfMeanMap_.find(filetag+"_ht") != sfMeanMap_.end() && !isSignal ) 
        {
            // Derive ht SF
            const double mean_ht = sfMeanMap_[filetag+"_ht"];
            htDerivedweight = (1/mean_ht)*htDerivedweightUncor;
            
            const double mean_ht_flat2000 = sfMeanMap_[filetag+"_ht_flat2000"];
            htDerivedweightFlat2000 = (1/mean_ht_flat2000)*htDerivedweightFlat2000Uncor;
            
            const double mean_ht_njet7 = sfMeanMap_[filetag+"_ht_njet7"];
            htDerivedweightNJet7 = (1/mean_ht_njet7)*htDerivedweightNJet7Uncor;
            
            const double mean_ht_MG = sfMeanMap_[filetag+"_ht_MG"];
            htDerivedweightMG = (1/mean_ht_MG)*htDerivedweightMGUncor;

            // Derive ht up and down variation on SF
            if( runYear == "2016" ) 
            {
                const double fit2NJetBin8 = 1.307*exp(-0.0003416*HT_trigger_pt30);
                const double fit2NJetBin567 = htScaleFactor(8, HT_trigger_pt30, runYear);
                const double ratioUp = fit2NJetBin8/fit2NJetBin567;
                const double ratioDown = fit2NJetBin567/fit2NJetBin8;
                
                const double fit2NJetBin8MG = 1.151*exp(-0.0001730*HT_trigger_pt30);
                const double fit2NJetBin567MG = htScaleFactor(8, HT_trigger_pt30, filetag);
                const double ratioUpMG = fit2NJetBin8/fit2NJetBin567;
                const double ratioDownMG = fit2NJetBin567/fit2NJetBin8;

                htScaleUp = htDerivedweight*ratioUp;
                htScaleDown = htDerivedweight*ratioDown;

                htScaleUpMG = htDerivedweightMG*ratioUpMG;
                htScaleDownMG = htDerivedweightMG*ratioDownMG;
            }
            else if( runYear == "2017" )
            {
                const double fit2NJetBin8 = 1.215*exp(-0.0002613*HT_trigger_pt30);
                const double fit2NJetBin567 = htScaleFactor(8, HT_trigger_pt30, runYear);
                const double ratioUp = fit2NJetBin8/fit2NJetBin567;
                const double ratioDown = fit2NJetBin567/fit2NJetBin8;
                
                const double fit2NJetBin8MG = 1.192*exp(-0.0002465*HT_trigger_pt30);
                const double fit2NJetBin567MG = htScaleFactor(8, HT_trigger_pt30, filetag);
                const double ratioUpMG = fit2NJetBin8/fit2NJetBin567;
                const double ratioDownMG = fit2NJetBin567/fit2NJetBin8;
                
                htScaleUp = htDerivedweight*ratioUp;
                htScaleDown = htDerivedweight*ratioDown;            
                
                htScaleUpMG = htDerivedweightMG*ratioUpMG;
                htScaleDownMG = htDerivedweightMG*ratioDownMG;
            }
        }

        tr.registerDerivedVar( "htDerivedweight"+myVarSuffix_, htDerivedweight);
        tr.registerDerivedVar( "htDerivedweightFlat2000"+myVarSuffix_, htDerivedweightFlat2000 );
        tr.registerDerivedVar( "htDerivedweightNJet7"+myVarSuffix_, htDerivedweightNJet7 );
        tr.registerDerivedVar( "htDerivedweightMG"+myVarSuffix_, htDerivedweightMG);

        tr.registerDerivedVar( "htDerivedweightUncor"+myVarSuffix_, htDerivedweightUncor);
        tr.registerDerivedVar( "htDerivedweightFlat2000Uncor"+myVarSuffix_, htDerivedweightFlat2000Uncor);
        tr.registerDerivedVar( "htDerivedweightNJet7Uncor"+myVarSuffix_, htDerivedweightNJet7Uncor );
        tr.registerDerivedVar( "htDerivedweightMGUncor"+myVarSuffix_, htDerivedweightMGUncor);

        tr.registerDerivedVar( "htScaleUp"+myVarSuffix_, htScaleUp);
        tr.registerDerivedVar( "htScaleDown"+myVarSuffix_, htScaleDown);
        tr.registerDerivedVar( "htScaleUpMG"+myVarSuffix_, htScaleUpMG);
        tr.registerDerivedVar( "htScaleDownMG"+myVarSuffix_, htScaleDownMG);

        //-----------------------------------------------------------------------------
        // Adding a scale factor for pileup 
        // For 2016: Grab the individual pileup weight from the histogram found in PileupHistograms_0121_69p2mb_pm4p6.root
        // For 2017: Grab the ratio from the histogram file and multiply this with the original weight 
        // ----------------------------------------------------------------------------
        
        const auto& puWeight   = tr.getVar<double>("puWeight");
        const auto& puSysUp    = tr.getVar<double>("puSysUp");
        const auto& puSysDown  = tr.getVar<double>("puSysDown");
        const auto& tru_npv    = tr.getVar<double>("TrueNumInteractions");
        double puWeightUnCorr  = 1.0, puSysUpUnCorr   = 1.0, puSysDownUnCorr = 1.0;

        if( runYear == "2016" || runYear == "2018" ) 
        {
            puWeightUnCorr = puSFHisto_->GetBinContent( findBin<TH1F, double>(puSFHisto_, tru_npv, "X") );
            puSysUpUnCorr = puSFUpHisto_->GetBinContent( findBin<TH1F, double>(puSFUpHisto_, tru_npv, "X") );
            puSysDownUnCorr = puSFDownHisto_->GetBinContent( findBin<TH1F, double>(puSFDownHisto_, tru_npv, "X") );
        }
        else if( runYear == "2017") 
        {
            puWeightUnCorr = puSFHisto_->GetBinContent( findBin<TH1F, double>(puSFHisto_, tru_npv, "X") )*puWeight; 
            puSysUpUnCorr = puSFUpHisto_->GetBinContent( findBin<TH1F, double>(puSFUpHisto_, tru_npv, "X") )*puWeight;
            puSysDownUnCorr = puSFDownHisto_->GetBinContent( findBin<TH1F, double>(puSFDownHisto_, tru_npv, "X") )*puWeight;
        }
        
        tr.registerDerivedVar( "puWeightUnCorr"+myVarSuffix_,  puWeightUnCorr);
        tr.registerDerivedVar( "puSysUpUnCorr"+myVarSuffix_,   puSysUpUnCorr);
        tr.registerDerivedVar( "puSysDownUnCorr"+myVarSuffix_, puSysDownUnCorr);

        // Adding correction to the pileup weight
        double puWeightCorr  = puWeight;
        double puSysUpCorr   = puSysUp;
        double puSysDownCorr = puSysDown;
        if( sfMeanMap_.find(filetag+"_pu") != sfMeanMap_.end() && !isSignal )
        {
            const double mean = sfMeanMap_[filetag+"_pu"];
            puWeightCorr     = (1/mean)*puWeightUnCorr;
            const double meanUp = sfMeanMap_[filetag+"_pu_Up"];
            puSysUpCorr      = (1/meanUp)*puSysUpUnCorr;
            const double meanDown = sfMeanMap_[filetag+"_pu_Down"];
            puSysDownCorr    = (1/meanDown)*puSysDownUnCorr;
        }

        tr.registerDerivedVar( "puWeightCorr"+myVarSuffix_, puWeightCorr);
        tr.registerDerivedVar( "puSysUpCorr"+myVarSuffix_, puSysUpCorr);
        tr.registerDerivedVar( "puSysDownCorr"+myVarSuffix_, puSysDownCorr);
        
        // --------------------------------------------------------------------------------------
        // Adding top pt reweighting for ttbar MC (Powheg)
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        // 13 TeV all combined
        // --------------------------------------------------------------------------------------
        double topPtScaleFactor = 1.0;
        auto* topPtVec = new std::vector<double>();
        if(filetag == "TT" || filetag == "2016_TT" || filetag.find("TTTo") != std::string::npos)
        {
            const double a=0.0615, b=-0.0005;
            auto SF = [&](const double pt){return exp(a + b*pt);};
            
            const auto& GenParticles        = tr.getVec<TLorentzVector>("GenParticles");
            const auto& GenParticles_PdgId  = tr.getVec<int>("GenParticles_PdgId");
            const auto& GenParticles_Status = tr.getVec<int>("GenParticles_Status"); 

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++)
            {
                if( abs(GenParticles_PdgId[gpi]) == 6 )
                {
                    topPtScaleFactor *= SF( GenParticles[gpi].Pt() );
                    topPtVec->push_back( GenParticles[gpi].Pt() );
                }
            }
            topPtScaleFactor = sqrt(topPtScaleFactor);
        }
        
        tr.registerDerivedVar( "topPtScaleFactor"+myVarSuffix_, topPtScaleFactor);
        tr.registerDerivedVec( "topPtVec"+myVarSuffix_, topPtVec);
        
        // --------------------------------------------------------------------------------------
        // Adding reweighting recipe to emulate Level 1 ECAL prefiring
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
        // --------------------------------------------------------------------------------------
        double prefiringScaleFactor = 1.0;
        if( runYear == "2017" )
        {
            const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);            
            const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);            
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodJets_pt30[ijet]) continue;
                const TLorentzVector& jet = Jets.at(ijet);
                const double weight = L1Prefireing_->GetBinContent(L1Prefireing_->FindBin(jet.Eta(), jet.Pt()));
                prefiringScaleFactor *= 1 - weight;
            }
        }
        tr.registerDerivedVar( "prefiringScaleFactor"+myVarSuffix_, prefiringScaleFactor);                    

        // Registering a variable that is the nominal total weight with lepton scale factor, btag scale factor, ht scale factor
        const auto& Weight         = tr.getVar<double>("Weight");
        const auto& bTagWeight     = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+myVarSuffix_);
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NGoodMuons     = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        
        double totalEventWeight    = -1.0;
        double totalEventWeightMG  = -1.0;
        if( NGoodElectrons == 1 ) 
        {
            totalEventWeight = Weight*bTagWeight*totGoodElectronSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
            totalEventWeightMG = Weight*bTagWeight*totGoodElectronSF*htDerivedweightMG*prefiringScaleFactor*puWeightCorr;
        }
        else if ( NGoodMuons == 1 ) 
        {
            totalEventWeight = Weight*bTagWeight*totGoodMuonSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
            totalEventWeightMG = Weight*bTagWeight*totGoodMuonSF*htDerivedweightMG*prefiringScaleFactor*puWeightCorr;
        }
        tr.registerDerivedVar( "totalEventWeight"+myVarSuffix_, totalEventWeight );
        tr.registerDerivedVar( "totalEventWeightMG"+myVarSuffix_, totalEventWeightMG );
    }

public:
    ScaleFactors( const std::string& runYear, const std::string& leptonFileName, const std::string& puFileName, const std::string& meanFileName, const std::string& myVarSuffix = "" )
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up ScaleFactors"<<std::endl;
        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?

        // Getting Lepton scale factor histograms
        TFile SFRootFile( leptonFileName.c_str() );
        TString eleSFHistoTightName, eleSFHistoIsoName, eleSFHistoRecoName, eleSFHistoTrigName;
        TString muSFHistoMediumName,  muSFHistoIsoName,  muSFHistoRecoName,  muSFHistoTrigName;
        if( runYear == "2016")
        {
            eleSFHistoTightName = "Run2016_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2016_Mini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2016_num_el_pt40_trig_5jCut_htCut_isoTrig";
            muSFHistoMediumName = "sf_mu_mediumID";
            muSFHistoIsoName = "sf_mu_mediumID_mini02";
            muSFHistoTrigName = "TrigEff_2016_num_mu_pt40_trig_5jCut_htCut_isoTrig";
        }
        else if( runYear == "2017")
        {
            eleSFHistoTightName = "Run2017_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2017_MVAVLooseTightIP2DMini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2017_num_el_pt40_trig_5jCut_htCut_isoTrig";
            muSFHistoMediumName = "NUM_MediumID_DEN_genTracks_pt_abseta";
            muSFHistoIsoName = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
            muSFHistoTrigName = "TrigEff_2017_num_mu_pt40_trig_5jCut_htCut_isoTrig";
        }
        else if( runYear == "2018")
        {
            eleSFHistoTightName = "Run2018_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2018_Mini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "";
            muSFHistoMediumName = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
            muSFHistoIsoName = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
            muSFHistoTrigName = "";
        }
        eleSFHistoTight_.reset( (TH2F*)SFRootFile.Get(eleSFHistoTightName) );
        eleSFHistoIso_.reset(   (TH2F*)SFRootFile.Get(eleSFHistoIsoName) );
        eleSFHistoReco_.reset(  (TH2F*)SFRootFile.Get(eleSFHistoRecoName) );
        eleSFHistoTrig_.reset(  (TH2F*)SFRootFile.Get(eleSFHistoTrigName) );
        muSFHistoMedium_.reset( (TH2F*)SFRootFile.Get(muSFHistoMediumName) );
        muSFHistoIso_.reset(    (TH2F*)SFRootFile.Get(muSFHistoIsoName) );
        muSFHistoTrig_.reset(   (TH2F*)SFRootFile.Get(muSFHistoTrigName) );
        if( runYear == "2016" ) 
        {
            eleSFHistoIP2D_.reset( (TH2F*)SFRootFile.Get("Run2016_MVAVLooseIP2D") );//In 2016, the isolation SF histogram is separate from the IP2D cut scale factor histogram.
            muSFHistoReco_.reset(  (TGraph*)SFRootFile.Get("ratio_eff_aeta_dr030e030_corr") ); //Only 2016 requires the track reconstruction efficiency.
        }
        SFRootFile.Close();

        // Getting mean of some scale factors to keep total number of events the same after apply the scale factor
        TFile SFMeanRootFile( meanFileName.c_str() );
        TIter next(SFMeanRootFile.GetListOfKeys());
        TKey* key;
        while(key = (TKey*)next())
        {
            std::shared_ptr<TH1> h( (TH1*)key->ReadObj() );
            std::string name( h->GetTitle() );
            sfMeanMap_.insert(std::pair<std::string, double>(name, h->GetMean()));
        }
        SFMeanRootFile.Close();

        // Getting the pile up histograms
        TFile puRootFile( puFileName.c_str() );
        if( runYear == "2016" || runYear == "2018") 
        {
            puSFHisto_.reset(     (TH1F*)puRootFile.Get("pu_weights_central") );
            puSFUpHisto_.reset(   (TH1F*)puRootFile.Get("pu_weights_up") );
            puSFDownHisto_.reset( (TH1F*)puRootFile.Get("pu_weights_down") );
        }
        else if( runYear == "2017")
        {
            puSFHisto_.reset(     (TH1F*)puRootFile.Get("pu_ratio_central") ); 
            puSFUpHisto_.reset(   (TH1F*)puRootFile.Get("pu_ratio_up") );
            puSFDownHisto_.reset( (TH1F*)puRootFile.Get("pu_ratio_down") ); 
        }
        puRootFile.Close();

        // Get the L1prefiring scale factor histogram
        TFile L1PrefiringFile("L1prefiring_jetpt_2017BtoF.root");
        L1Prefireing_.reset( (TH2F*)L1PrefiringFile.Get("L1prefiring_jetpt_2017BtoF") );
        L1PrefiringFile.Close();
    }

    ~ScaleFactors() 
    {
    }

    void operator()(NTupleReader& tr)
    {
        scaleFactors(tr);
    }
};

#endif
