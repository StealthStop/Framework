#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
class ScaleFactors
{
private:

    TH2F *eleSFHistoTight_;        
    TH2F *eleSFHistoIso_;
    TH2F *eleSFHistoReco_;
    TH2F *muSFHisto_;
    TGraph *muSFHistoReco_;

    void scaleFactors(NTupleReader& tr)
    {

        //First for PDF Uncertainties
        const auto& scaleWeights    = tr.getVec<double>("ScaleWeights");
        const auto& PDFweights      = tr.getVec<double>("PDFweights");

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
        if( !std::isfinite( scaleWeightNominal ) ) {
            scaleWeightNominal = 1.0;
            scaleWeightUpperBound = 1.0;
            scaleWeightLowerBound = 1.0;
        }

        tr.registerDerivedVar("scaleWeightUp",      scaleWeightUpperBound);
        tr.registerDerivedVar("scaleWeightDown",    scaleWeightLowerBound);
        tr.registerDerivedVar("scaleWeightNom",     scaleWeightNominal);

        //Now calculate the PDF scale factor and uncertainty based on the 100 different replica values stored in PDFweights using envelope method and the median

        auto PDFWeightMax            = std::max_element( std::begin(PDFweights), std::end(PDFweights) );
        auto PDFWeightMin            = std::min_element( std::begin(PDFweights), std::end(PDFweights) );

        double PDFWeightUpperBound   = *PDFWeightMax;
        double PDFWeightLowerBound   = *PDFWeightMin;

        const double reqCL           = 0.68; //Choose a confidence level for the uncertainty
        std::vector<double> sortedPDFWeights = PDFweights; //Cannot sort a constant
        std::sort( sortedPDFWeights.begin() + 1, sortedPDFWeights.end() );
        
        const int upper = std::round( 0.5 * (1 + reqCL) * 100.0 );
        const int lower = 1 + std::round( 0.5 * (1 - reqCL) * 100.0 );

        double central  = 0.5*( sortedPDFWeights[50] + sortedPDFWeights[51] ); //Exactly 100 entries
        double errminus = central - sortedPDFWeights[lower];
        double errplus  = sortedPDFWeights[upper] - central;
        double errsymm  = 0.5*( errplus + errminus );

        double NNPDF_from_median_up = central + errplus;
        NNPDF_from_median_up = ( ( NNPDF_from_median_up/central ) > 2.0 ) ?  1.0 : ( ( ( NNPDF_from_median_up/central ) < -2.0 ) ? 1.0 : NNPDF_from_median_up/central );
        
        double NNPDF_from_median_down = central - errplus;
        NNPDF_from_median_down = NNPDF_from_median_down/central > 2.0 ? 1.0 : NNPDF_from_median_down/central < -2.0 ? 1.0 : NNPDF_from_median_down/central;
        
        if( !std::isfinite(central) ) {
            NNPDF_from_median_up    = 1.0;
            NNPDF_from_median_down  = 1.0;
            central                 = 1.0;
        }

        tr.registerDerivedVar( "PDFweightUp",   NNPDF_from_median_up );
        tr.registerDerivedVar( "PDFweightDown", NNPDF_from_median_down );
        tr.registerDerivedVar( "PDFweightNom",  central );
        
        //Adding code for implementing electron scale factors


//        TFile SFRootFile( SFRootFileName_.c_str() );

        const auto& electrons           = tr.getVec<TLorentzVector>("Electrons");
        const auto& goodElectrons       = tr.getVec<bool>("GoodElectrons");
        
        double totGoodElectronSF        = 1.0;
        double totGoodElectronSFPErr2   = 0.0;
        double totGoodElectronSFErr     = 0.0;
        double totGoodElectronSF_Up     = 1.0;
        double totGoodElectronSF_Down   = 1.0;
        
        //Loop through the electrons
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) {

            //If it is not a good lepton, give scale factor of 1.0
            if( !goodElectrons.at(iel) ) continue;

            else {
                //Get the scale factor from the rootfile
                double elpt     = electrons.at(iel).Pt();
                double eleta    = electrons.at(iel).Eta();
                int xbin = 0, ybin = 0, ptbin = 0, etabin = 0;

                for( unsigned int ixbin = 0; ixbin < eleSFHistoTight_->GetNbinsX()+1; ixbin++ ) {
                    double tempxBinEdgeMax = (double) eleSFHistoTight_->GetXaxis()->GetBinUpEdge(ixbin);
                    if( elpt < tempxBinEdgeMax )  {
                        xbin = ixbin; 
                        break;
                    }
                }

                for( unsigned int iybin = 0; iybin < eleSFHistoTight_->GetNbinsY()+1; iybin++ ) {
                    double tempyBinEdgeMax = (double)eleSFHistoTight_->GetYaxis()->GetBinUpEdge(iybin);
                    if( std::fabs(eleta) < tempyBinEdgeMax ) {
                        ybin = iybin;
                        break;
                    }
                }
                
                if( xbin == 0 && elpt > 200.0) xbin = eleSFHistoTight_->GetNbinsX(); //std::cout<<elpt<<std::endl;
                if( ybin == 0 ) std::cerr<<"Invalid eta stored for a good electron!"<<std::endl;

                for( unsigned int iptbin = 0; iptbin < eleSFHistoReco_->GetNbinsY()+1; iptbin++ ) {
                    double tempPtBinEdgeMax = (double) eleSFHistoReco_->GetYaxis()->GetBinUpEdge(iptbin);
                    if( elpt < tempPtBinEdgeMax ) {
                        ptbin = iptbin;
                        break;
                    }
                }

                for( unsigned int ietabin = 0; ietabin < eleSFHistoReco_->GetNbinsX()+1; ietabin++ ) {
                    double tempEtaBinEdgeMax = (double) eleSFHistoReco_->GetXaxis()->GetBinUpEdge(ietabin);
                    if( eleta < tempEtaBinEdgeMax ) {
                        etabin = ietabin;
                        break;
                    }
                }
                
                if( ptbin == 0 && elpt > 500.0) ptbin = eleSFHistoReco_->GetNbinsY(); //std::cout<<elpt<<std::endl;
                if( etabin == 0 ) std::cerr<<"Invalid eta stored for a good electron!"<<std::endl;


                if( xbin != 0 && ybin != 0 && ptbin != 0 && etabin != 0 ) {
                    double eleTightSF       = eleSFHistoTight_->GetBinContent( xbin, ybin );
                    double eleTightSFErr    = eleSFHistoTight_->GetBinError( xbin, ybin );
                    double eleTightPErr     = eleTightSFErr/eleTightSF;

                    double eleIsoSF         = eleSFHistoIso_->GetBinContent( xbin, ybin );
                    double eleIsoSFErr      = eleSFHistoIso_->GetBinError( xbin, ybin );
                    double eleIsoPErr       = eleIsoSFErr/eleIsoSF;

                    double eleRecoSF        = eleSFHistoReco_->GetBinContent( etabin, ptbin );
                    double eleRecoSFErr     = eleSFHistoReco_->GetBinError( etabin, ptbin );
                    double eleRecoPErr      = eleRecoSFErr/eleRecoSF;

                    //The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    double eleTotSF         = eleTightSF * eleIsoSF * eleRecoSF; 
                    double eleTotPErr       = std::sqrt( eleTightPErr*eleTightPErr + eleIsoPErr*eleIsoPErr + eleRecoPErr*eleRecoPErr );
                    double eleTotSFErr      = eleTotPErr * eleTotSF;

                    //if( eleTotSF < 0.1 ) {
                    //    std::cout<<"EL pt: "<<elpt<<"; eta:"<<eleta<<"; "<<xbin<<" "<<ybin<<" "<<etabin<<" "<<ptbin<<" "<<eleSFHistoTight->GetBinContent(xbin,ybin)<<" "<<eleSFHistoIso->GetBinContent(xbin,ybin)<<" "<<eleSFHistoReco->GetBinContent(etabin, ptbin)<<std::endl;
                    //}
                    
                    totGoodElectronSF       *= eleTotSF; 
                    totGoodElectronSFPErr2  += eleTotPErr * eleTotPErr;
                }
            }
        }

        totGoodElectronSFErr    = std::sqrt(totGoodElectronSFPErr2) * totGoodElectronSF;
        totGoodElectronSF_Up    = totGoodElectronSF + totGoodElectronSFErr;
        totGoodElectronSF_Down  = totGoodElectronSF - totGoodElectronSFErr;

        tr.registerDerivedVar( "totGoodElectronSF",         totGoodElectronSF );
        tr.registerDerivedVar( "totGoodElectronSFErr",      totGoodElectronSFErr );
        tr.registerDerivedVar( "totGoodElectronSF_Up",      totGoodElectronSF_Up );
        tr.registerDerivedVar( "totGoodElectronSF_Down",    totGoodElectronSF_Down );

        //Adding code for implementing muon scale factors

        const auto& muons           = tr.getVec<TLorentzVector>("Muons");
        const auto& goodMuons       = tr.getVec<bool>("GoodMuons");

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
                
                int xbin = 0, ybin = 0, etabin = 0;

                for( unsigned int ixbin = 0; ixbin < muSFHisto_->GetNbinsX()+1; ixbin++ ) {
                    double tempxBinEdgeMax = muSFHisto_->GetXaxis()->GetBinUpEdge(ixbin);
                    if( mupt < tempxBinEdgeMax )  {
                        xbin = ixbin; 
                        break;
                    }
                }

                for( unsigned int iybin = 0; iybin < muSFHisto_->GetNbinsY()+1; iybin++ ) {
                    double tempyBinEdgeMax = muSFHisto_->GetYaxis()->GetBinUpEdge(iybin);
                    if( mueta < tempyBinEdgeMax ) {
                        ybin = iybin;
                        break;
                    }
                }
               
                if( xbin == 0 && mupt > 200.0) xbin = muSFHisto_->GetNbinsX(); //If the muon does not have a pT smaller than the max bin edge, it must have a pT greater than 200, which according to the Twiki, means we use the largest pT scale factor (until further notice).
                if( ybin == 0 ) std::cerr<<"Invalid eta stored for a good muon!"<<std::endl;//Our good leptons should not have an absolute value of eta greater than that of the max bin of the histogram
                
                //std::cout<<"MU pt: "<<mupt<<" ; eta: "<<mueta<<"; "<<xbin<<" "<<ybin<<" "<<muSFHisto->GetBinContent(xbin,ybin)<<" "<<muSFHistoReco->Eval(mueta)<<";"<<muSFHisto->GetBinContent(xbin,ybin)*muSFHistoReco->Eval(mueta)<<std::endl;

                if( xbin != 0 && ybin != 0) {

                    //The SUSLepton Twiki claims that the errors in the histogrm are purely statistical and can be ignored and recommends a 3% error for each leg
                    double muIDSF           = muSFHisto_->GetBinContent( xbin, ybin );
                    double muIDSFPErr       = .03;
                    double muIDSFErr        = muIDSF * muIDSFPErr;
                    
                    //For the general track reconstruction they claim that the systematics for the systematic still need to be finalized.
                    double muRecoSF         = muSFHistoReco_->Eval( mueta );

                    double muTotSF          = muIDSF * muRecoSF;
                    double muTotSFPErr      = muIDSFPErr; //This is a place holder for if and when we get the reconstruction scale factor systematic

                    totGoodMuonSF           *= muTotSF;
                    totGoodMuonSFPErr2      += muTotSFPErr * muTotSFPErr;
                }
            }
        }

        totGoodMuonSFErr    = std::sqrt(totGoodMuonSFPErr2) * totGoodMuonSF;
        totGoodMuonSF_Up    = totGoodMuonSF + totGoodMuonSFErr;
        totGoodMuonSF_Down  = totGoodMuonSF - totGoodMuonSFErr;

        tr.registerDerivedVar( "totGoodMuonSF",         totGoodMuonSF );
        tr.registerDerivedVar( "totGoodMuonSFErr",      totGoodMuonSFErr );
        tr.registerDerivedVar( "totGoodMuonSF_Up",      totGoodMuonSF_Up );
        tr.registerDerivedVar( "totGoodMuonSF_Down",    totGoodMuonSF_Down );
    }

public:
    ScaleFactors( std::string SFRootFileName = "2016ScaleFactorHistos.root" )
        : eleSFHistoTight_(nullptr)
        , eleSFHistoIso_(nullptr)
        , eleSFHistoReco_(nullptr)
        , muSFHisto_(nullptr)
        , muSFHistoReco_(nullptr)
          
    {
        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?
        TFile SFRootFile( SFRootFileName.c_str() );
        
        eleSFHistoTight_       = (TH2F*)SFRootFile.Get("GsfElectronToCutBasedSpring15T");
        eleSFHistoIso_         = (TH2F*)SFRootFile.Get("MVAVLooseElectronToMini");
        eleSFHistoReco_        = (TH2F*)SFRootFile.Get("EGamma_SF2D");
        
        muSFHisto_             = (TH2F*)SFRootFile.Get("sf_mu_mediumID_mini02");
        muSFHistoReco_         = (TGraph*)SFRootFile.Get("ratio_eff_aeta_dr030e030_corr");

        SFRootFile.Close();
    }

    ~ScaleFactors() {
    }

    void operator()(NTupleReader& tr)
    {
        scaleFactors(tr);
    }
};

#endif
