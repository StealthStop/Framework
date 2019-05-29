#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
#include "TFile.h"
#include "TKey.h"

class ScaleFactors
{
private:
    std::string myVarSuffix_;

    TH2F* eleSFHistoTight_;        
    TH2F* eleSFHistoIso_;
    TH2F* eleSFHistoIP2D_;
    TH2F* eleSFHistoReco_;
    TH2F* eleSFHistoTrig_;
    TH2F* muSFHistoMedium_;
    TH2F* muSFHistoIso_;
    TH2F* muSFHistoTrig_;
    TGraph* muSFHistoReco_;
    TH1F* puSFHisto_;
    TH1F* puSFUpHisto_;
    TH1F* puSFDownHisto_;
    std::shared_ptr<TH2F> L1Prefireing_;
    std::map<std::string, double> sfMeanMap_;

    void scaleFactors(NTupleReader& tr)
    {
        // --------------------------------------------------------------------------------------
        // First for PDF Uncertainties
        // --------------------------------------------------------------------------------------
        const auto& scaleWeights    = tr.getVec<double>("ScaleWeights");
        const auto& PDFweights      = tr.getVec<double>("PDFweights");
        const auto& filetag         = tr.getVar<std::string>("filetag");

        //Following the example in SusyAnaTools PDFUncertainty.h, the scale weights are calculated using the envelope method and we ignore all anti-correlated variations (5 and 7)

        std::vector<double>  myScaleWeights;
        /*std::cout<<"Length of scale weight vector (must be 9): "<<scaleWeights.size()<<std::endl;
        for( unsigned int i = 0; i < scaleWeights.size(); i++ ) {
            std::cout<<"Scale weight vector entry number "<<i+1<<": "<<scaleWeights.at(i)<<std::endl;
        }*/
        if( scaleWeights.size() == 9 ) {//If there are not exactly 9 scale factors, then the vector  was filled incorrectly
            myScaleWeights.push_back( scaleWeights.at(1) );
            myScaleWeights.push_back( scaleWeights.at(2) );
            myScaleWeights.push_back( scaleWeights.at(3) );
            myScaleWeights.push_back( scaleWeights.at(4) );
            myScaleWeights.push_back( scaleWeights.at(6) );
            myScaleWeights.push_back( scaleWeights.at(8) );
        }
        else {
            myScaleWeights.clear();
            myScaleWeights.resize(6, 1.0);
        }

        auto scaleWeightMax          = std::max_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        auto scaleWeightMin          = std::min_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        double scaleWeightUpperBound = *scaleWeightMax;
        double scaleWeightLowerBound = *scaleWeightMin;
        double scaleWeightNominal    = scaleWeights.size() == 9 ?  scaleWeights.at(0) : 1.0;

        //TODO: There are some NaN/Inf values in the Diboson channel - still need to figure out why this is an issue.
        //std::cout<<scaleWeightUpperBound<<" "<<scaleWeightNominal<<" "<<scaleWeightLowerBound<<std::endl;
        if( !std::isfinite(scaleWeightNominal) || !std::isfinite(scaleWeightUpperBound) || !std::isfinite(scaleWeightLowerBound) ) 
        {
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
        // Store different partonshower uncertainties in more easy to parse way
        // Note: not all samples have these weights stored, give them default value of 1. 
        // --------------------------------------------------------------------------------------        
        double PSweight_ISRUp = 1.; 
        double PSweight_ISRDown = 1.; 
        double PSweight_FSRUp = 1.; 
        double PSweight_FSRDown = 1.; 
        double PSweight_ISRUp_2 = 1.; 
        double PSweight_ISRDown_2 = 1.; 
        double PSweight_FSRUp_2 = 1.; 
        double PSweight_FSRDown_2 = 1.; 
        if(tr.hasVar("PSweights")){
            const auto& PSweights         = tr.getVec<double>("PSweights");
            if(&PSweights != nullptr && PSweights.size() >= 12) // should have size of 14, but just put 12 or more to be able to use the sample with the bug
            {
                // Get nominal one so we can normalize it
                double MEweight = PSweights.at(0);
                // reduced variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/sqrt(2) and sqrt(2)
                PSweight_ISRUp = PSweights.at(2)/MEweight;
                PSweight_FSRUp = (PSweights.at(3)/MEweight < 10.0) ? PSweights.at(3)/MEweight : 1.0;
                PSweight_ISRDown = PSweights.at(4)/MEweight;
                PSweight_FSRDown = (PSweights.at(5)/MEweight < 10.0) ? PSweights.at(5)/MEweight : 1.0;
                // nominal variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/2 and 2
                PSweight_ISRUp_2 = PSweights.at(6)/MEweight;
                PSweight_FSRUp_2 = (PSweights.at(7)/MEweight < 10.0) ? PSweights.at(7)/MEweight : 1.0;
                PSweight_ISRDown_2 = PSweights.at(8)/MEweight;
                PSweight_FSRDown_2 = (PSweights.at(9)/MEweight < 10.0) ? PSweights.at(9)/MEweight : 1.0;
            }
            //else
            //{
            //    std::cout << "PS weights not in expected format, " << PSweights.size() << " weights found (12 or 14 expected). " << std::endl;
            //}
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
        // Now calculate the PDF scale factor and uncertainty based on the 100 different replica values stored in PDFweights using envelope method and the median
        // --------------------------------------------------------------------------------------
        double central = 1.0;
        double NNPDF_from_median_up = 1.0;
        double NNPDF_from_median_down = 1.0;
        if(PDFweights.size() > 0)
        {
            auto PDFWeightMax            = std::max_element( std::begin(PDFweights), std::end(PDFweights) );
            auto PDFWeightMin            = std::min_element( std::begin(PDFweights), std::end(PDFweights) );

            double PDFWeightUpperBound   = *PDFWeightMax;
            double PDFWeightLowerBound   = *PDFWeightMin;

            const double reqCL           = 0.68; //Choose a confidence level for the uncertainty
            std::vector<double> sortedPDFWeights = PDFweights; //Cannot sort a constant
            std::sort( sortedPDFWeights.begin() + 1, sortedPDFWeights.end() );
        
            const int upper = std::round( 0.5 * (1 + reqCL) * 100.0 );
            const int lower = 1 + std::round( 0.5 * (1 - reqCL) * 100.0 );

            central  = 0.5*( sortedPDFWeights[50] + sortedPDFWeights[51] ); //Exactly 100 entries
            double errminus = central - sortedPDFWeights[lower];
            double errplus  = sortedPDFWeights[upper] - central;
            double errsymm  = 0.5*( errplus + errminus );

            NNPDF_from_median_up = central + errplus;
            NNPDF_from_median_up = ( ( NNPDF_from_median_up/central ) > 2.0 ) ?  1.0 : ( ( ( NNPDF_from_median_up/central ) < -2.0 ) ? 1.0 : NNPDF_from_median_up/central );
        
            NNPDF_from_median_down = central - errplus;
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
        // Adding code for implementing electron scale factors
        // --------------------------------------------------------------------------------------
        const auto& electrons           = tr.getVec<TLorentzVector>("Electrons");
        const auto& goodElectrons       = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        
        double totGoodElectronSF        = 1.0;
        double totGoodElectronSFPErr2   = 0.0;
        double totGoodElectronSFErr     = 0.0;
        double totGoodElectronSF_Up     = 1.0;
        double totGoodElectronSF_Down   = 1.0;
        
        //Loop through the electrons
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) {

            //If it is not a good lepton, give scale factor of 1.0
            if( !goodElectrons.at(iel) ) continue;

            else {//Not sure if this is necessary, but do not want to make too many changes
                //Get the scale factor from the rootfile
                double elpt     = electrons.at(iel).Pt();
                double eleta    = electrons.at(iel).Eta();
                int xbinElTight = 0, ybinElTight = 0, xbinElIso = 0, ybinElIso = 0, xbinElReco = 0, ybinElReco = 0, xbinElTrig = 0, ybinElTrig = 0;

                if( filetag.find("2017") != std::string::npos ) { //If this is the year 2017
                    
                    //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for Tight ID ( has 180 bins total in 2D parameter space )
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoTight_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoTight_->GetXaxis()->GetBinUpEdge( ixbin );
                        if( eleta < tempxBinEdgeMax ) {
                            xbinElTight = ixbin;
                            break;
                        }
                    }
                    for( unsigned int iybin = 0; iybin < eleSFHistoTight_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoTight_->GetYaxis()->GetBinUpEdge ( iybin );
                        if( elpt < tempyBinEdgeMax ) {
                            ybinElTight = iybin;
                            break;
                        }
                    }
                    
                    if( xbinElTight == 0 ) std::cerr<<"El Tight Histo: Invalid eta stored for a good electron!"<<std::endl;
                    //If the pt of the lepton is larger than the maximum value of the scale factor chart, default to the scale factor in the largest bin.
                    if( ybinElTight == 0 && elpt > 500.0) ybinElTight = eleSFHistoTight_->GetNbinsY(); //std::cout<<elpt<<std::endl;

                    //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for MiniIso of 0.1  ( has 210 bins total in 2D parameter space )
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoIso_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoIso_->GetXaxis()->GetBinUpEdge( ixbin );
                        if( eleta < tempxBinEdgeMax ) {
                            xbinElIso = ixbin;
                            break;
                        }
                    }

                    for( unsigned int iybin = 0; iybin < eleSFHistoIso_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoIso_->GetYaxis()->GetBinUpEdge ( iybin );
                        if( elpt < tempyBinEdgeMax ) {
                            ybinElIso = iybin;
                            break;
                        }
                    }
                    if( xbinElIso == 0 ) std::cerr<<"El Iso Histo: Invalid eta stored for a good electron!"<<std::endl;
                    //If the pt of the lepton is larger than the maximum value of the scale factor chart, default to the scale factor in the largest bin.
                    if( ybinElIso == 0 && elpt > 500.0) ybinElIso = eleSFHistoIso_->GetNbinsY(); //std::cout<<elpt<<std::endl;
                    
                    //Find the bin indices (binned by x: eta and y: pt ) for the 2017 scale factor for Data/MC comparison (reco eff) ( has 144 bins total in 2D parameter space )
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoReco_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoReco_->GetXaxis()->GetBinUpEdge( ixbin );
                        if( eleta < tempxBinEdgeMax ) {
                            xbinElReco = ixbin;
                            break;
                        }
                    }

                    for( unsigned int iybin = 0; iybin < eleSFHistoReco_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoReco_->GetYaxis()->GetBinUpEdge ( iybin );
                        if( elpt < tempyBinEdgeMax ) {
                            ybinElReco = iybin;
                            break;
                        }
                    }
                    if( xbinElReco == 0 ) std::cerr<<"El Reco Histo: Invalid eta stored for a good electron!"<<std::endl;
                    //If the pt of the lepton is larger than the maximum value of the scale factor chart, default to the scale factor in the largest bin.
                    if( ybinElReco == 0 && elpt > 500.0) ybinElReco = eleSFHistoReco_->GetNbinsY(); //std::cout<<elpt<<std::endl;
                    
                    //Find the bin indices (binned by x: pt and y: eta ) for the 2017 scale factor for Data/MC trig efficiency
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoTrig_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoTrig_->GetXaxis()->GetBinUpEdge( ixbin );
                        if( elpt < tempxBinEdgeMax ) {
                            xbinElTrig = ixbin;
                            break;
                        }
                    }

                    for( unsigned int iybin = 0; iybin < eleSFHistoTrig_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoTrig_->GetYaxis()->GetBinUpEdge ( iybin );
                        if( eleta < tempyBinEdgeMax ) {
                            ybinElTrig = iybin;
                            break;
                        }
                    }
                    
                    if( xbinElTrig == 0 && elpt > 200.0) xbinElTrig = eleSFHistoTrig_->GetNbinsX(); //std::cout<<elpt<<std::endl;
                    if( ybinElTrig == 0 ) std::cerr<<"El Trig Histo: Invalid eta stored for a good electron!"<<std::endl;
                    //If the pt of the lepton is larger than the maximum value of the scale factor chart, default to the scale factor in the largest bin.
                        
                    //std::cout<<"EL pt: "<<elpt<<"; eta:"<<eleta<<"; "<<xbinElTight<<" "<<ybinElTight<<" "<<xbinElIso<<" "<<ybinElIso<<" "<<xbinElReco<<" "<<ybinElReco<<" "<<eleSFHistoTight_->GetBinContent(xbinElTight,ybinElTight)<<" "<<eleSFHistoIso_->GetBinContent(xbinElIso,ybinElIso)<<" "<<eleSFHistoReco_->GetBinContent(xbinElReco, ybinElReco)<<std::endl;
                }//END OF 2017 Loop

                else { //For 2016 electrons
                
                    //Find the bin indices (binned by x: eta and y: pt ) for the 2016 scale factor for Tight ID ( has 30 bins total in 2D parameter space )
                    //This is the same index as for the Iso ID
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoTight_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoTight_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( eleta < tempxBinEdgeMax )  {
                            xbinElTight = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < eleSFHistoTight_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double)eleSFHistoTight_->GetYaxis()->GetBinUpEdge(iybin);
                        if( elpt < tempyBinEdgeMax ) {
                            ybinElTight = iybin;
                            break;
                        }
                    }
                    
                    if( xbinElTight == 0 ) std::cerr<<"Invalid eta stored for a good electron!"<<std::endl;
                    if( ybinElTight == 0 && elpt > 500.0) ybinElTight = eleSFHistoTight_->GetNbinsY(); 
                    //Since the binning for MiniIso < 0.1 is the same as that for Tight ID, we will use the same values (values initialized for uniformity later on). Same for the IP2D (hence you do not need an extra set of variables).
                    xbinElIso = xbinElTight;
                    ybinElIso = ybinElTight;
    
                    //Find the bin indices (binned by x: eta and y: pt ) for the 2016 scale factor for Data/MC differences (reco eff) ( has 90 bins total in 2D parameter space )
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoReco_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoReco_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( eleta < tempxBinEdgeMax ) {
                            xbinElReco = ixbin;
                            break;
                        }
                    }
                    
                    for( unsigned int iybin = 0; iybin < eleSFHistoReco_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoReco_->GetYaxis()->GetBinUpEdge(iybin);
                        if( elpt < tempyBinEdgeMax ) {
                            ybinElReco = iybin;
                            break;
                        }
                    }
    
                    if( xbinElReco == 0 ) std::cerr<<"El Reco Histo: Invalid eta stored for a good electron!"<<std::endl;
                    if( ybinElReco == 0 && elpt > 500.0) ybinElReco = eleSFHistoReco_->GetNbinsY(); //std::cout<<elpt<<std::endl;
                    
                    //Find the bin indices (binned by x: pt and y: eta ) for the 2016 scale factor for Data/MC trig efficiency
                    for( unsigned int ixbin = 0; ixbin < eleSFHistoTrig_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = (double) eleSFHistoTrig_->GetXaxis()->GetBinUpEdge( ixbin );
                        if( elpt < tempxBinEdgeMax ) {
                            xbinElTrig = ixbin;
                            break;
                        }
                    }

                    for( unsigned int iybin = 0; iybin < eleSFHistoTrig_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = (double) eleSFHistoTrig_->GetYaxis()->GetBinUpEdge ( iybin );
                        if( eleta < tempyBinEdgeMax ) {
                            ybinElTrig = iybin;
                            break;
                        }
                    }
                    
                    if( xbinElTrig == 0 && elpt > 200.0) xbinElTrig = eleSFHistoTrig_->GetNbinsX(); //std::cout<<elpt<<std::endl;
                    if( ybinElTrig == 0 ) std::cerr<<"El Trig Histo: Invalid eta stored for a good electron!"<<std::endl;
                } //END of 2016 Loop

                if( xbinElTight != 0 && ybinElTight != 0 && xbinElIso != 0 && ybinElIso != 0 && xbinElReco != 0 && ybinElReco != 0 ) {

                    double eleTightSF       = eleSFHistoTight_->GetBinContent( xbinElTight, ybinElTight );
                    double eleTightSFErr    = eleSFHistoTight_->GetBinError( xbinElTight, ybinElTight );
                    double eleTightPErr     = eleTightSFErr/eleTightSF;

                    double eleIsoSF         = eleSFHistoIso_->GetBinContent( xbinElIso, ybinElIso );
                    double eleIsoSFErr      = eleSFHistoIso_->GetBinError( xbinElIso, ybinElIso );
                    double eleIsoPErr       = eleIsoSFErr/eleIsoSF;

                    double eleRecoSF        = eleSFHistoReco_->GetBinContent( xbinElReco, ybinElReco );
                    double eleRecoSFErr     = eleSFHistoReco_->GetBinError( xbinElReco, ybinElReco );
                    double eleRecoPErr      = eleRecoSFErr/eleRecoSF;
                    
                    double eleTrigSF        = eleSFHistoTrig_->GetBinContent( xbinElTrig, ybinElTrig );
                    double eleTrigSFErr     = eleSFHistoTrig_->GetBinError( xbinElTrig, ybinElTrig );
                    double eleTrigPErr      = eleTrigSFErr/eleTrigSF;

                    //The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    if( filetag.find("2016") != std::string::npos ) { //If this is the year 2016, we need to add the IP2D histogram scale factors into the Iso scale factor
                        double eleIP2DSF    = eleSFHistoIP2D_->GetBinContent( xbinElIso, ybinElIso );
                        double eleIP2DSFErr = eleSFHistoIP2D_->GetBinError( xbinElIso, ybinElIso );
                        double eleIP2DPErr  = eleIP2DSFErr/eleIP2DSF;

                        eleIsoSF            = eleIsoSF*eleIP2DSF;
                        eleIsoPErr          = std::sqrt( eleIsoPErr*eleIsoPErr + eleIP2DPErr*eleIP2DPErr );
                        eleIsoSFErr         = eleIsoPErr*eleIsoSF;
                    }
                    
                    double eleTotSF         = eleTightSF * eleIsoSF * eleRecoSF * eleTrigSF; 
                    double eleTotPErr       = std::sqrt( eleTightPErr*eleTightPErr + eleIsoPErr*eleIsoPErr + eleRecoPErr*eleRecoPErr + eleTrigPErr*eleTrigPErr);
                    double eleTotSFErr      = eleTotPErr * eleTotSF;

                    //if( eleTotSF < 0.1 ) {
                    //    std::cout<<"EL pt: "<<elpt<<"; eta:"<<eleta<<"; "<<xbinElTight<<" "<<ybinElTight<<" "<<xbinElIso<<" "<<ybinElIso<<" "<<xbinElReco<<" "<<ybinElReco<<" "<<eleSFHistoTight_->GetBinContent(xbinElTight,ybinElTight)<<" "<<eleSFHistoIso_->GetBinContent(xbinElIso,ybinElIso)<<" "<<eleSFHistoReco_->GetBinContent(xbinElReco, ybinElReco)<<std::endl;
                    //}
                    
                    totGoodElectronSF       *= eleTotSF; 
                    totGoodElectronSFPErr2  += eleTotPErr * eleTotPErr;
                }
            }
        }

        totGoodElectronSFErr    = std::sqrt(totGoodElectronSFPErr2) * totGoodElectronSF;
        totGoodElectronSF_Up    = totGoodElectronSF + totGoodElectronSFErr;
        totGoodElectronSF_Down  = totGoodElectronSF - totGoodElectronSFErr;

        tr.registerDerivedVar( "totGoodElectronSF"+myVarSuffix_,         totGoodElectronSF );
        tr.registerDerivedVar( "totGoodElectronSFErr"+myVarSuffix_,      totGoodElectronSFErr );
        tr.registerDerivedVar( "totGoodElectronSF_Up"+myVarSuffix_,      totGoodElectronSF_Up );
        tr.registerDerivedVar( "totGoodElectronSF_Down"+myVarSuffix_,    totGoodElectronSF_Down );

        // --------------------------------------------------------------------------------------
        // Adding code for implementing muon scale factors
        // --------------------------------------------------------------------------------------
        const auto& muons           = tr.getVec<TLorentzVector>("Muons");
        const auto& goodMuons       = tr.getVec<bool>("GoodMuons"+myVarSuffix_);

        double totGoodMuonSF        = 1.0;
        double totGoodMuonSFPErr2   = 0.0;
        double totGoodMuonSFErr     = 0.0;
        double totGoodMuonSF_Up     = 1.0;
        double totGoodMuonSF_Down   = 1.0;

        //Loop through the muons
        for( unsigned int imu = 0; imu < muons.size(); imu++ ) {
            
            //If it is not a good lepton, no need to calculate scale factor since we are not using it in the analysis
            if( !goodMuons.at(imu) ) continue;

            else {
                //Get the scale factor from the rootfile
                double mupt = muons.at(imu).Pt();
                double mueta = std::fabs( muons.at(imu).Eta() );
                
                int xbinMuMedium = 0, ybinMuMedium = 0, xbinMuIso = 0, ybinMuIso = 0, xbinMuTrig = 0, ybinMuTrig = 0, etabin = 0;
 
                if( filetag.find("2017") != std::string::npos ) { //If this is the year 2017
                    //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2017 scale factor for Medium ID (has only 24 bins)
                    for( unsigned int ixbin = 0; ixbin < muSFHistoMedium_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = muSFHistoMedium_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( mupt < tempxBinEdgeMax )  {
                            xbinMuMedium = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < muSFHistoMedium_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = muSFHistoMedium_->GetYaxis()->GetBinUpEdge(iybin);
                        if( std::fabs( mueta ) < tempyBinEdgeMax ) { //Histogram is binned by absolute value of eta!
                            ybinMuMedium = iybin;
                            break;
                        }
                    }
                    if( xbinMuMedium == 0 && mupt > 120.0) xbinMuMedium = muSFHistoMedium_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                    if( ybinMuMedium == 0 ) std::cerr<<"Mu ID Histo: Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram
                    //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2017 scale factor for MiniIso of 0.2 (has only 20 bins)
                    for( unsigned int ixbin = 0; ixbin < muSFHistoIso_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = muSFHistoIso_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( mupt < tempxBinEdgeMax )  {
                            xbinMuIso = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < muSFHistoIso_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = muSFHistoIso_->GetYaxis()->GetBinUpEdge(iybin);
                        if( std::fabs( mueta ) < tempyBinEdgeMax ) { //Histogram is binned by absolute value of eta!
                            ybinMuIso = iybin;
                            break;
                        }
                    }
                    if( xbinMuIso == 0 && mupt > 120.0) xbinMuIso = muSFHistoIso_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                    if( ybinMuIso == 0 ) std::cerr<<"Mu Iso Histo: Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram
                    //Find the bin indices (binned by x: pt and y: eta ) for the 2017 scale factor for Data/MC Trigger Efficiency
                    for( unsigned int ixbin = 0; ixbin < muSFHistoTrig_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = muSFHistoTrig_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( mupt < tempxBinEdgeMax )  {
                            xbinMuTrig = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < muSFHistoTrig_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = muSFHistoTrig_->GetYaxis()->GetBinUpEdge(iybin);
                        if( mueta < tempyBinEdgeMax ) {
                            ybinMuTrig = iybin;
                            break;
                        }
                    }
                    if( xbinMuTrig == 0 && mupt > 200.0) xbinMuTrig = muSFHistoTrig_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                    if( ybinMuTrig == 0 ) std::cerr<<"Mu Trig Histo: Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram
                }//END of 2017 loop

                else { //If the year is 2016 loop
                    
                    //Find the bin indices (binned by x: pt and y: abs(eta) ) for the 2016 scale factor for both Medium ID and MiniIso of 0.2 (has only 20 bins) - same binning!
                    for( unsigned int ixbin = 0; ixbin < muSFHistoMedium_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = muSFHistoMedium_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( mupt < tempxBinEdgeMax )  {
                            xbinMuMedium = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < muSFHistoMedium_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = muSFHistoMedium_->GetYaxis()->GetBinUpEdge(iybin);
                        if( std::fabs( mueta ) < tempyBinEdgeMax ) {
                            ybinMuMedium = iybin;
                            break;
                        }
                    }
                   
                    if( xbinMuMedium == 0 && mupt > 200.0) xbinMuMedium = muSFHistoMedium_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                    if( ybinMuMedium == 0 ) std::cerr<<"Mu Iso Histo: Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram

                    //For consistency with 2017, copy the values over to these variables.
                    xbinMuIso = xbinMuMedium;
                    ybinMuIso = ybinMuMedium;

                    //Find the bin indices (binned by x: pt and y: eta ) for the 2016 scale factor for Data/MC Trigger Efficiency
                    for( unsigned int ixbin = 0; ixbin < muSFHistoTrig_->GetNbinsX()+1; ixbin++ ) {
                        double tempxBinEdgeMax = muSFHistoTrig_->GetXaxis()->GetBinUpEdge(ixbin);
                        if( mupt < tempxBinEdgeMax )  {
                            xbinMuTrig = ixbin; 
                            break;
                        }
                    }
    
                    for( unsigned int iybin = 0; iybin < muSFHistoTrig_->GetNbinsY()+1; iybin++ ) {
                        double tempyBinEdgeMax = muSFHistoTrig_->GetYaxis()->GetBinUpEdge(iybin);
                        if( mueta < tempyBinEdgeMax ) {
                            ybinMuTrig = iybin;
                            break;
                        }
                    }
                    if( xbinMuTrig == 0 && mupt > 200.0) xbinMuTrig = muSFHistoTrig_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                    if( ybinMuTrig == 0 ) std::cerr<<"Mu Trig Histo: Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram
                }//END of 2016 loop
                //std::cout<<"MU pt: "<<mupt<<" ; eta: "<<mueta<<"; "<<xbin<<" "<<ybin<<" "<<muSFHisto->GetBinContent(xbin,ybin)<<" "<<muSFHistoReco->Eval(mueta)<<";"<<muSFHisto->GetBinContent(xbin,ybin)*muSFHistoReco->Eval(mueta)<<std::endl;

                if( xbinMuMedium != 0 && ybinMuMedium != 0 && xbinMuIso != 0 && ybinMuIso != 0 ) {

                    //The SUSLepton Twiki claims that the errors in the histogrm are purely statistical and can be ignored and recommends a 3% error for each leg (ID+IP+ISO)
                    double muMediumSF           = muSFHistoMedium_->GetBinContent( xbinMuMedium, ybinMuMedium );
                    double muIsoSF              = muSFHistoIso_->GetBinContent( xbinMuIso, ybinMuIso );
                    double muTrigSF             = muSFHistoTrig_->GetBinContent( xbinMuTrig, ybinMuTrig );
                    double muTrigSFErr          = muSFHistoTrig_->GetBinError( xbinMuTrig, ybinMuTrig );
                    double muTrigSFPErr         = muTrigSFErr/muTrigSF;

                    double muTotSF              = muMediumSF * muIsoSF * muTrigSF;
                    double muTotSFPErr2         = 0.03*0.03 + muTrigSFPErr*muTrigSFPErr;

                    if( filetag.find("2017") == std::string::npos ) {
                    //For the general track reconstruction they claim that the errors for the systematic still need to be finalized - does not seem to have been finalized as of Dec 2018
                    //This reconstruction value only exists for 2016 - SUS SF people say the 3% will include the reco scale factor uncertainty for now
                        double muRecoSF         = muSFHistoReco_->Eval( mueta );
                        muTotSF                 = muTotSF * muRecoSF;
                    }
                    totGoodMuonSF           *= muTotSF;
                    totGoodMuonSFPErr2      += muTotSFPErr2;
                }
            }
        }

        totGoodMuonSFErr    = std::sqrt(totGoodMuonSFPErr2) * totGoodMuonSF;
        totGoodMuonSF_Up    = totGoodMuonSF + totGoodMuonSFErr;
        totGoodMuonSF_Down  = totGoodMuonSF - totGoodMuonSFErr;

        tr.registerDerivedVar( "totGoodMuonSF"+myVarSuffix_,         totGoodMuonSF );
        tr.registerDerivedVar( "totGoodMuonSFErr"+myVarSuffix_,      totGoodMuonSFErr );
        tr.registerDerivedVar( "totGoodMuonSF_Up"+myVarSuffix_,      totGoodMuonSF_Up );
        tr.registerDerivedVar( "totGoodMuonSF_Down"+myVarSuffix_,    totGoodMuonSF_Down );

        // --------------------------------------------------------------------------------------
        // Adding a scale factor that corrects the disagreement between data and MC for Ht
        // --------------------------------------------------------------------------------------
        // Moved the filetag since I need it earlier also for the lepton scale factors
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& isSignal = tr.getVar<bool>("isSignal");

        auto htScaleFactor = [](int nJets, double HT, std::string filetag) 
        { 
            double norm     = 0.0;
            double expo     = 0.0;
            if( filetag.find("2016") != std::string::npos ) { //If it is for 2016 MC sample
                //No PU version
                //norm = 0.05669*nJets + 0.8391;
                //expo = (-0.04318*nJets - 0.03314)/1000;
                
                //PU version
                //norm = 0.03175*nJets + 0.9504;
                //expo = (-0.02100*nJets - 0.1031 )/1000;
                
                //New PU version
                norm = 0.03422*nJets + 0.9367;
                expo = (-0.02310*nJets - 0.0940 )/1000;
            }
            else {
                //No PU version
                //norm = 0.0149284*nJets+1.00437;
                //expo = (-0.0019273*nJets-0.134854)/1000;
                
                //PU version
                //norm = 0.03897*nJets + 0.8993;
                //expo = (-0.02851*nJets - 0.04683 )/1000;
                
                //New PU version
                norm = 0.02565*nJets + 0.9635;
                expo = (-0.01418*nJets - 0.1101 )/1000;
            }
            return norm*exp( expo*HT ); 
        };
        auto htScaleFactorFlat2000 = [](int nJets, double HT, std::string filetag)
        {
            double norm     = 0.0;
            double expo     = 0.0;
            
            if( filetag.find("2016") != std::string::npos ) { //If it is for 2016 MC sample
                //No PU version
                //norm = 0.05669*nJets + 0.8391;
                //expo = (-0.04318*nJets - 0.03314)/1000;
                
                //PU version
                //norm = 0.03175*nJets + 0.9504;
                //expo = (-0.02100*nJets - 0.1031 )/1000;
                
                //New PU version
                norm = 0.03422*nJets + 0.9367;
                expo = (-0.02310*nJets - 0.0940 )/1000;
                
            }
            else {
                //No PU version
                //norm = 0.0149284*nJets+1.00437;
                //expo = (-0.0019273*nJets-0.134854)/1000;
                
                //PU version
                //norm = 0.03897*nJets + 0.8993;
                //expo = (-0.02851*nJets - 0.04683 )/1000;
                
                //New PU version
                norm = 0.02565*nJets + 0.9635;
                expo = (-0.01418*nJets - 0.1101 )/1000;
            }
            
            if( HT > 2000 ) {
                return norm*exp( expo*2000.00 );
            }
            
            else {
                return norm*exp( expo*HT ); 
            }
        };
        auto htScaleFactorNJet7 = [](double HT, std::string filetag)
        {
            double norm     = 0.0;
            double expo     = 0.0;
            
            if( filetag.find("2016") != std::string::npos ) { //If it is for 2016 MC sample
                //No PU version
                //norm = 0.05669*7 + 0.8391;
                //expo = (-0.04318*7 - 0.03314)/1000;
                
                //PU version
                //norm = 0.03175*7 + 0.9504;
                //expo = (-0.02100*7 - 0.1031 )/1000;
                
                //New PU version
                norm = 0.03422*7 + 0.9367;
                expo = (-0.02310*7 - 0.0940 )/1000;
            }
            else {
                //No PU version
                //norm = 0.0149284*7+1.00437;
                //expo = (-0.0019273*7-0.134854)/1000;
                
                //PU version
                //norm = 0.03897*7 + 0.8993;
                //expo = (-0.02851*7 - 0.04683 )/1000;
                
                //New PU version
                norm = 0.02565*7 + 0.9635;
                expo = (-0.01418*7 - 0.1101 )/1000;
            }
            
            return norm*exp( expo*HT ); 
        };

        double htDerivedweight = 1.0;
        double htDerivedweightFlat2000 = 1.0;
        double htDerivedweightNJet7 = 1.0;

        double htDerivedweightUncor = htScaleFactor(NGoodJets_pt30, HT_trigger_pt30, filetag);
        double htDerivedweightFlat2000Uncor = htScaleFactorFlat2000(NGoodJets_pt30, HT_trigger_pt30, filetag);
        double htDerivedweightNJet7Uncor = htScaleFactorNJet7(HT_trigger_pt30, filetag);

        double htScaleUp = 1.0;
        double htScaleDown = 1.0;

        if( sfMeanMap_.find(filetag+"_ht") != sfMeanMap_.end() && !isSignal ) 
        {
            // Derive ht SF
            const double mean_ht = sfMeanMap_[filetag+"_ht"];
            htDerivedweight = (1/mean_ht)*htDerivedweightUncor;
            
            const double mean_ht_flat2000 = sfMeanMap_[filetag+"_ht_flat2000"];
            htDerivedweightFlat2000 = (1/mean_ht_flat2000)*htDerivedweightFlat2000Uncor;
            
            const double mean_ht_njet7 = sfMeanMap_[filetag+"_ht_njet7"];
            htDerivedweightNJet7 = (1/mean_ht_njet7)*htDerivedweightNJet7Uncor;

            // Derive ht up and down variation on SF
            if( filetag.find("2016") != std::string::npos ) {
                const double fit2NJetBin8 = 1.307*exp(-0.0003416*HT_trigger_pt30);
                const double fit2NJetBin567 = htScaleFactor(8, HT_trigger_pt30, filetag);
                const double ratioUp = fit2NJetBin8/fit2NJetBin567;
                const double ratioDown = fit2NJetBin567/fit2NJetBin8;

                htScaleUp = htDerivedweight*ratioUp;
                htScaleDown = htDerivedweight*ratioDown;
            }
            else {
                const double fit2NJetBin8 = 1.215*exp(-0.0002613*HT_trigger_pt30);
                const double fit2NJetBin567 = htScaleFactor(8, HT_trigger_pt30, filetag);
                const double ratioUp = fit2NJetBin8/fit2NJetBin567;
                const double ratioDown = fit2NJetBin567/fit2NJetBin8;
                
                htScaleUp = htDerivedweight*ratioUp;
                htScaleDown = htDerivedweight*ratioDown;
            
            }
        }

        tr.registerDerivedVar( "htDerivedweight"+myVarSuffix_, htDerivedweight);
        tr.registerDerivedVar( "htDerivedweightFlat2000"+myVarSuffix_, htDerivedweightFlat2000 );
        tr.registerDerivedVar( "htDerivedweightNJet7"+myVarSuffix_, htDerivedweightNJet7 );

        tr.registerDerivedVar( "htDerivedweightUncor"+myVarSuffix_, htDerivedweightUncor);
        tr.registerDerivedVar( "htDerivedweightFlat2000Uncor"+myVarSuffix_, htDerivedweightFlat2000Uncor);
        tr.registerDerivedVar( "htDerivedweightNJet7Uncor"+myVarSuffix_, htDerivedweightNJet7Uncor );

        tr.registerDerivedVar( "htScaleUp"+myVarSuffix_, htScaleUp);
        tr.registerDerivedVar( "htScaleDown"+myVarSuffix_, htScaleDown);


        //-----------------------------------------------------------------------------
        //
        // For 2016: Grab the individual pileup weight from the histogram found in PileupHistograms_0121_69p2mb_pm4p6.root
        // For 2017: Grab the ratio from the histogram file and multiply this with the original weight 
        //
        // ----------------------------------------------------------------------------
        
        const auto puWeight         = tr.getVar<double>("puWeight");
        const auto puSysUp          = tr.getVar<double>("puSysUp");
        const auto puSysDown        = tr.getVar<double>("puSysDown");

        const auto& tru_npv          = tr.getVar<double>("TrueNumInteractions");
        double puWeightUnCorr        = 1.0;
        double puSysUpUnCorr         = 1.0;
        double puSysDownUnCorr       = 1.0;

        if( filetag.find("2016") != std::string::npos ) { //If this is the year 2016
            if( tru_npv < puSFHisto_->GetBinLowEdge( puSFHisto_->GetNbinsX()+1 ) ) {
                puWeightUnCorr = puSFHisto_->GetBinContent( puSFHisto_->GetXaxis()->FindBin(tru_npv) );
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFHisto_"<<std::endl;
                puWeightUnCorr = puSFHisto_->GetBinContent( puSFHisto_->GetNbinsX() );
            }
            
            if( tru_npv < puSFUpHisto_->GetBinLowEdge( puSFUpHisto_->GetNbinsX()+1 ) ) {
                puSysUpUnCorr = puSFUpHisto_->GetBinContent( puSFUpHisto_->GetXaxis()->FindBin(tru_npv) );
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFUpHisto_"<<std::endl;
                puSysUpUnCorr = puSFUpHisto_->GetBinContent( puSFUpHisto_->GetNbinsX() );
            }
            
            if( tru_npv < puSFDownHisto_->GetBinLowEdge( puSFDownHisto_->GetNbinsX()+1 ) ) {
                puSysDownUnCorr = puSFDownHisto_->GetBinContent( puSFDownHisto_->GetXaxis()->FindBin(tru_npv) );
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFDownHisto_"<<std::endl;
                puSysDownUnCorr = puSFDownHisto_->GetBinContent( puSFDownHisto_->GetNbinsX() );
            }
        }

        else { //If this is the year 2017
            if( tru_npv< puSFHisto_->GetBinLowEdge( puSFHisto_->GetNbinsX() + 1 ) ) {
                puWeightUnCorr = puSFHisto_->GetBinContent( puSFHisto_->GetXaxis()->FindBin(tru_npv) )*puWeight; 
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFHisto_"<<std::endl;
                puWeightUnCorr = puSFHisto_->GetBinContent( puSFHisto_->GetNbinsX() )*puWeight;
            }
            if( tru_npv< puSFUpHisto_->GetBinLowEdge( puSFUpHisto_->GetNbinsX() + 1 ) ) {
                puSysUpUnCorr = puSFUpHisto_->GetBinContent( puSFUpHisto_->GetXaxis()->FindBin(tru_npv) )*puWeight; 
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFUpHisto_"<<std::endl;
                puSysUpUnCorr = puSFUpHisto_->GetBinContent( puSFUpHisto_->GetNbinsX() )*puSysUp;
            }
            if( tru_npv< puSFDownHisto_->GetBinLowEdge( puSFDownHisto_->GetNbinsX() + 1 ) ) {
                puSysDownUnCorr = puSFDownHisto_->GetBinContent( puSFDownHisto_->GetXaxis()->FindBin(tru_npv) )*puWeight; 
            }
            else {
                std::cerr<<"The true num of interactions is larger than the maximum number of bins in puSFDownHisto_"<<std::endl;
                puSysDownUnCorr = puSFDownHisto_->GetBinContent( puSFDownHisto_->GetNbinsX() )*puSysDown;
            }
        }
        
        tr.registerDerivedVar( "puWeightUnCorr"+myVarSuffix_,  puWeightUnCorr);
        tr.registerDerivedVar( "puSysUpUnCorr"+myVarSuffix_,   puSysUpUnCorr);
        tr.registerDerivedVar( "puSysDownUnCorr"+myVarSuffix_, puSysDownUnCorr);
        
        // --------------------------------------------------------------------------
        //
        // Adding correction to the pileup weight
        //
        //--------------------------------------------------------------------------

        double puWeightCorr         = puWeight;
        double puSysUpCorr          = puSysUp;
        double puSysDownCorr        = puSysDown;

        if( sfMeanMap_.find(filetag+"_pu") != sfMeanMap_.end() && !isSignal )
        {
            const double mean = sfMeanMap_[filetag+"_pu"];
            puWeightCorr            = (1/mean)*puWeightUnCorr;
            const double meanUp = sfMeanMap_[filetag+"_pu_Up"];
            puSysUpCorr             = (1/meanUp)*puSysUpUnCorr;
            const double meanDown = sfMeanMap_[filetag+"_pu_Down"];
            puSysDownCorr           = (1/meanDown)*puSysDownUnCorr;
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
        if(filetag.find("2017") != std::string::npos)
        {
            const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);            
            const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);            
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodJets_pt30[ijet]) continue;
                TLorentzVector jet = Jets.at(ijet);
                const auto& eta = jet.Eta();
                const auto& pt = jet.Pt();
                double weight = L1Prefireing_->GetBinContent(L1Prefireing_->FindBin(eta, pt));
                prefiringScaleFactor *= 1 - weight;
            }
        }
        tr.registerDerivedVar( "prefiringScaleFactor"+myVarSuffix_, prefiringScaleFactor);                    
        // Registering a variable that is the nominal total weight with lepton scale factor, btag scale factor, ht scale factor
        const auto& Weight              = tr.getVar<double>("Weight");
        const auto& bTagWeight          = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+myVarSuffix_);
        int    NGoodElectrons           = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        int    NGoodMuons               = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        
        double totalEventWeight         = -1.0;

        if( NGoodElectrons == 1 ) {
            totalEventWeight = Weight*bTagWeight*totGoodElectronSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
        }
        else if ( NGoodMuons == 1 ) {
            totalEventWeight = Weight*bTagWeight*totGoodMuonSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
        }

        tr.registerDerivedVar( "totalEventWeight"+myVarSuffix_, totalEventWeight );
    }
    
public:
    ScaleFactors( const std::string& SFRootFileName = "2016ScaleFactorHistos.root", const std::string& puRootFileName = "PileupHistograms_0121_69p2mb_pm4p6.root", const std::string& SFMeanRootFileName = "allInOne_SFMean.root", const std::string& myVarSuffix = "" )
        : myVarSuffix_(myVarSuffix)
        , eleSFHistoTight_(nullptr)
        , eleSFHistoIso_(nullptr)
        , eleSFHistoIP2D_(nullptr)
        , eleSFHistoReco_(nullptr)
        , eleSFHistoTrig_(nullptr)
        , muSFHistoMedium_(nullptr)
        , muSFHistoIso_(nullptr)
        , muSFHistoReco_(nullptr)
        , muSFHistoTrig_(nullptr)
        , puSFHisto_(nullptr)
        , puSFUpHisto_(nullptr)
        , puSFDownHisto_(nullptr)
    {
        std::cout<<"Setting up ScaleFactors"<<std::endl;
        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?
        TFile SFRootFile( SFRootFileName.c_str() );

        TString eleSFHistoTightName = ( SFRootFileName.find("2017") != std::string::npos ) ? "Run2017_CutBasedTightNoIso94XV2" : "Run2016_CutBasedTightNoIso94XV2";
        TString eleSFHistoIsoName = ( SFRootFileName.find("2017") != std::string::npos ) ? "Run2017_MVAVLooseTightIP2DMini" : "Run2016_Mini";
        TString eleSFHistoTrigName = ( SFRootFileName.find("2017") != std::string::npos ) ? "TrigEff_2017_num_el_pt40_trig_5jCut_htCut_isoTrig" : "TrigEff_2016_num_el_pt40_trig_5jCut_htCut_isoTrig";
        TString muSFHistoMediumName = ( SFRootFileName.find("2017") != std::string::npos ) ? "NUM_MediumID_DEN_genTracks_pt_abseta" : "sf_mu_mediumID"; 
        TString muSFHistoIsoName = ( SFRootFileName.find("2017") != std::string::npos ) ? "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta" : "sf_mu_mediumID_mini02";
        TString muSFHistoTrigName = ( SFRootFileName.find("2017") != std::string::npos ) ? "TrigEff_2017_num_mu_pt40_trig_5jCut_htCut_isoTrig" : "TrigEff_2016_num_mu_pt40_trig_5jCut_htCut_isoTrig";
        
        eleSFHistoTight_        = (TH2F*)SFRootFile.Get(eleSFHistoTightName);
        eleSFHistoIso_          = (TH2F*)SFRootFile.Get(eleSFHistoIsoName);
        eleSFHistoReco_         = (TH2F*)SFRootFile.Get("EGamma_SF2D"); //The name for this histogram is the same for both 2016 and 2017
        eleSFHistoTrig_         = (TH2F*)SFRootFile.Get(eleSFHistoTrigName);
        
        muSFHistoMedium_        = (TH2F*)SFRootFile.Get(muSFHistoMediumName);
        muSFHistoIso_           = (TH2F*)SFRootFile.Get(muSFHistoIsoName);
        muSFHistoTrig_          = (TH2F*)SFRootFile.Get(muSFHistoTrigName);
        if ( SFRootFileName.find("2017") == std::string::npos ) {
            muSFHistoReco_         = (TGraph*)SFRootFile.Get("ratio_eff_aeta_dr030e030_corr"); //Only 2016 requires the track reconstruction efficiency.
            eleSFHistoIP2D_        = (TH2F*)SFRootFile.Get("Run2016_MVAVLooseIP2D");//In 2016, the isolation SF histogram is separate from the IP2D cut scale factor histogram.
        }

        SFRootFile.Close();

        TFile SFMeanRootFile( SFMeanRootFileName.c_str() );
        TIter next(SFMeanRootFile.GetListOfKeys());
        TKey* key;
        while(key = (TKey*)next())
        {
            std::shared_ptr<TH1> h( (TH1*)key->ReadObj() );
            std::string name( h->GetTitle() );
            sfMeanMap_.insert(std::pair<std::string, double>(name, h->GetMean()));
        }
        SFMeanRootFile.Close();

        TFile puRootFile( puRootFileName.c_str() );
        
        if( puRootFileName.find("PileupHistograms_0121_69p2mb_pm4p6.root") != std::string::npos ) {
            puSFHisto_          = (TH1F*)puRootFile.Get("pu_weights_central");
            puSFUpHisto_        = (TH1F*)puRootFile.Get("pu_weights_up");
            puSFDownHisto_      = (TH1F*)puRootFile.Get("pu_weights_down");
        }
        else {
            puSFHisto_          = (TH1F*)puRootFile.Get("pu_ratio_central"); 
            puSFUpHisto_        = (TH1F*)puRootFile.Get("pu_ratio_up");
            puSFDownHisto_      = (TH1F*)puRootFile.Get("pu_ratio_down"); 
        }

        puRootFile.Close();

        TFile L1PrefiringFile("L1prefiring_jetpt_2017BtoF.root");
        L1Prefireing_.reset( (TH2F*)L1PrefiringFile.Get("L1prefiring_jetpt_2017BtoF") );
        L1PrefiringFile.Close();
    }

    ~ScaleFactors() {
    }

    void operator()(NTupleReader& tr)
    {
        scaleFactors(tr);
    }
};

#endif
