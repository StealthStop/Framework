#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
class ScaleFactors
{
private:

    std::string eleSFRootFileName_;
    std::string muSFRootFileName_;

    void scaleFactors(NTupleReader& tr)
    {
        // Get needed branches

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
        const auto& electrons           = tr.getVec<TLorentzVector>("Electrons");
        const auto& goodElectrons       = tr.getVec<bool>("GoodElectrons");
        
        auto* electronsSF               = new std::vector<double>;
        auto* electronsSFErr            = new std::vector<double>;
        auto* electronsSF_Up            = new std::vector<double>;
        auto* electronsSF_Down          = new std::vector<double>;
        
        double totGoodElectronSF        = 1.0;
        double totGoodElectronSFErr     = 1.0;
        double totGoodElectronSF_Up     = 1.0;
        double totGoodElectronSF_Down   = 1.0;

        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?

        TFile eleSFRootFile( muSFRootFileName_.c_str() );
        
        TH2F *eleSFHistoTight = (TH2F*)eleSFRootFile.Get("GsfElectronToCutBasedSpring15T");
        TH2F *eleSFHistoIso   = (TH2F*)eleSFRootFile.Get("MVAVLooseElectronToMini");
        TH2F *eleSFHistoReco  = (TH2F*)eleSFRootFile.Get("EGamma_SF2D");

        //Loop through the electrons
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) {
            //If it is not a good lepton, give scale factor of 1.0
            if( !goodElectrons.at(iel) ) electronSFs->push_back(1.0);

            else {
                //Get the scale factor from the rootfile
                double elpt     = electrons.at(iel).Pt();
                double eleta    = electrons.at(iel).Eta();
                int xbin = 0, ybin = 0, ptbin = 0, etabin = 0;

                for( unsigned int ixbin = 0; ixbin < eleSFHistoTight->GetNbinsX()+1; ixbin++ ) {
                    double tempxBinEdgeMax = (double) eleSFHistoTight->GetXaxis()->GetBinUpEdge(ixbin);
                    if( elpt < tempxBinEdgeMax )  {
                        xbin = ixbin; 
                        break;
                    }
                }

                for( unsigned int iybin = 0; iybin < eleSFHistoTight->GetNbinsY()+1; iybin++ ) {
                    double tempyBinEdgeMax = (double)eleSFHistoTight->GetYaxis()->GetBinUpEdge(iybin);
                    if( std::fabs(eleta) < tempyBinEdgeMax ) {
                        ybin = iybin;
                        break;
                    }
                }
                
                if( xbin == 0 && elpt > 200.0) xbin = eleSFHistoTight->GetNbinsX(); //std::cout<<elpt<<std::endl;
                if( ybin == 0 ) std::cerr<<"Invalid eta stored for a good electron!"<<std::endl;

                for( unsigned int iptbin = 0; iptbin < eleSFHistoReco->GetNbinsY()+1; iptbin++ ) {
                    double tempPtBinEdgeMax = (double) eleSFHistoReco->GetYaxis()->GetBinUpEdge(iptbin);
                    if( elpt < tempPtBinEdgeMax ) {
                        ptbin = iptbin;
                        break;
                    }
                }

                for( unsigned int ietabin = 0; ietabin < eleSFHistoReco->GetNbinsX()+1; ietabin++ ) {
                    double tempEtaBinEdgeMax = (double) eleSFHistoReco->GetXaxis()->GetBinUpEdge(ietabin);
                    if( eleta < tempEtaBinEdgeMax ) {
                        etabin = ietabin;
                        break;
                    }
                }
                
                if( ptbin == 0 && elpt > 500.0) ptbin = eleSFHistoReco->GetNbinsY(); //std::cout<<elpt<<std::endl;
                if( etabin == 0 ) std::cerr<<"Invalid eta stored for a good electron!"<<std::endl;


                if( xbin != 0 && ybin != 0 && ptbin != 0 && etabin != 0 ) {
                    double eleTightSF       = eleSFHistoTight->GetBinContent( xbin, ybin );
                    double eleTightSFErr    = eleSFHistoTight->GetBinError( xbin, ybin );
                    double eleTightPErr     = eleTightSFErr/eleTightSF;

                    double eleIsoSF         = eleSFHistoIso->GetBinContent( xbin, ybin );
                    double eleIsoSFErr      = eleSFHistoIso->GetBinError( xbin, ybin );
                    double eleIsoPErr       = eleIsoSFErr/eleIsoSF;

                    double eleRecoSF        = eleSFHistoReco->GetBinContent( etabin, ptbin );
                    double eleRecoSFErr     = eleSFHistoReco->GetBinError( etabin, ptbin );
                    double eleRecoPErr      = eleRecoSFErr/eleRecoSF;

                    //The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    double eleTotSF         = eleTightSF * eleIsoSF * eleRecoSF; 
                    double eleTotPErr       = std::sqrt( eleTightPErr*eleTightPErr + eleIsoPErr*eleIsoPErr + eleRecoPErr*eleRecoPErr );
                    double eleTotSFErr      = eleTotPErr * eleTotSF;

                    double eleTotSF_Up      = eleTotSF + eleTotSFErr;
                    double eleTotSF_Down    = eleTotSF - eleTotSFErr;

                    if( eleTotSF < 0.1 ) {
                        std::cout<<"EL pt: "<<elpt<<"; eta:"<<eleta<<"; "<<xbin<<" "<<ybin<<" "<<etabin<<" "<<ptbin<<" "<<eleSFHistoTight->GetBinContent(xbin,ybin)<<" "<<eleSFHistoIso->GetBinContent(xbin,ybin)<<" "<<eleSFHistoReco->GetBinContent(etabin, ptbin)<<std::endl;
                    }
                    
                    electronsSF->push_back(eleTotSF);
                    electronsSFErr->push_back(eleTotSFErr);
                    electronsSF_Up->push_back(eleTotSF_Up);
                    electronsSF_Down->push_back(eleTotSF_Down);
                    
                }
            }
        }

        if( electronsSF.size() != 0 ) {
            for( unsigned int iSF = 0; iSF < electronsSF.size(); iSF++ + {
                
            }
        }

        tr.registerDerivedVar( "totGoodElectronSF",         totGoodElectronSF );
        tr.registerDerivedVar( "totGoodElectronSFErr",      totGoodElectronSFErr );
        tr.registerDerivedVar( "totGoodElectronSF_Up",      totGoodElectronSF_Up );
        tr.registerDerivedVar( "totGoodElectronSF_Down",    totGoodElectronSF_Down );

        tr.registerDerivedVec( "electronsSF",       electronsSF );
        tr.registerDerivedVec( "electronsSFErr",    electronsSFErr );
        tr.registerDerivedVec( "electronsSF_Up",    electronsSF_Up );
        tr.registerDerivedVec( "electronsSF_Down",  electronsSF_Down );

        //Adding code for implementing muon scale factors
        TFile muSFRootFile( muSFRootFileName_.c_str() );

        TH2F *muSFHisto     = (TH2F*)muSFRootFile.Get("sf_mu_mediumID_mini02");
        TGraph *muSFHistoReco = (TGraph*)muSFRootFile.Get("ratio_eff_aeta_dr030e030_corr");

        const auto& muons           = tr.getVec<TLorentzVector>("Muons");
        const auto& goodMuons       = tr.getVec<bool>("GoodMuons");

        auto* muonSFs               = new std::vector<double>;
        double leadingMuonSF = -1.0;

        //Loop through the muons
        for( unsigned int imu = 0; imu < muons.size(); imu++ ) {
            //If it is not a good lepton, give scale factor of 1.0
            if( !goodMuons.at(imu) ) muonSFs->push_back(1.0);

            else {
                //Get the scale factor from the rootfile
                double mupt = muons.at(imu).Pt();
                double mueta = std::fabs( muons.at(imu).Eta() );
                
                int xbin = 0, ybin = 0, etabin = 0;

                for( unsigned int ixbin = 0; ixbin < muSFHisto->GetNbinsX()+1; ixbin++ ) {
                    double tempxBinEdgeMax = muSFHisto->GetXaxis()->GetBinUpEdge(ixbin);
                    if( mupt < tempxBinEdgeMax )  {
                        xbin = ixbin; 
                        break;
                    }
                }

                for( unsigned int iybin = 0; iybin < muSFHisto->GetNbinsY()+1; iybin++ ) {
                    double tempyBinEdgeMax = muSFHisto->GetYaxis()->GetBinUpEdge(iybin);
                    if( mueta < tempyBinEdgeMax ) {
                        ybin = iybin;
                        break;
                    }
                }
                
                
                if( xbin == 0 && mupt > 200.0) xbin = muSFHisto->GetNbinsX(); //std::cout<<elpt<<std::endl;
                if( ybin == 0 ) std::cerr<<"Invalid eta stored for a good muon!"<<std::endl;
                
                //std::cout<<"MU pt: "<<mupt<<" ; eta: "<<mueta<<"; "<<xbin<<" "<<ybin<<" "<<muSFHisto->GetBinContent(xbin,ybin)<<" "<<muSFHistoReco->Eval(mueta)<<";"<<muSFHisto->GetBinContent(xbin,ybin)*muSFHistoReco->Eval(mueta)<<std::endl;

                if( xbin != 0 && ybin != 0) {

                    double muTotSF = muSFHisto->GetBinContent( xbin, ybin ) * muSFHistoReco->Eval( mueta );
                    std::cout<<muTotSF<<std::endl;
                    muonSFs->push_back(muTotSF);
                    
                    //If leadingElectronSF is -1.0, then set it to the first electron SF obtained (ranked by pT)
                    if( leadingMuonSF < 0.0) {
                        leadingMuonSF = muTotSF;
                    }
                }
            }
        }
        
        tr.registerDerivedVar( "leadGoodMuonSF", leadingMuonSF );
        tr.registerDerivedVec( "muonSFs", muonSFs );
    }

public:
    ScaleFactors( std::string eleSFRootFileName = "ElectronScaleFactors_Run2016.root", std::string muSFRootFileName = "allINone_leptonSF_Moriond17.root")
        : eleSFRootFileName_(eleSFRootFileName)
        , muSFRootFileName_(muSFRootFileName)
    {}

    void operator()(NTupleReader& tr)
    {
        scaleFactors(tr);
    }
};

#endif
