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
    bool firstPrint_;
    bool printMeanError_;
    bool printGetBinError_;
    std::set<std::string> binError_message_;

    std::shared_ptr<TH2F> eleSFHistoTight_;        
    std::shared_ptr<TH2F> eleSFHistoIso_;
    std::shared_ptr<TH2F> eleSFHistoIP2D_;
    std::shared_ptr<TH2F> eleSFHistoReco_;
    std::shared_ptr<TH2F> eleSFHistoTrig_;
    std::shared_ptr<TH2F> muSFHistoMedium_;
    std::shared_ptr<TH2F> muSFHistoIso_;
    std::shared_ptr<TH2F> muSFHistoTrig_;
    std::shared_ptr<TH2F> nimuSFHistoTrig_;
    std::shared_ptr<TGraph> muSFHistoReco_;
    std::shared_ptr<TH2F> L1Prefireing_;
    std::map<std::string, double> sfMeanMap_;
    std::shared_ptr<TH2F> jetSFHistoTrigName_2bCut_;
    std::shared_ptr<TH2F> jetSFHistoTrigName_3bCut_;
    std::shared_ptr<TH2F> jetSFHistoTrigName_ge4bCut_;

    template<typename T> std::shared_ptr<T>& getHisto(TFile& f, std::shared_ptr<T>& h, const TString& name)
    {        
        if(name != "") h.reset( static_cast<T*>(f.Get(name)) );   
        else std::cerr<<utility::color("Warning: A needed scale factor histogram is set to nullptr and using 1.0 as the default", "red")<<std::endl;
        return h;
    }

    template<typename Th, typename Tb> int findBin(const std::shared_ptr<Th>& h, const Tb v, const std::string& axis, const std::string& message = "")
    {
        int bin = -1;
        if(h)
        {
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
        }
        if(bin == -1 && firstPrint_)
        {
            binError_message_.insert(message); 
            printGetBinError_ = true;
        }
        return bin;
    }

    double htScaleFactor(const int nJets, const double HT, const std::string& runYear) const
    {
        // All values updated for v1.0 on October 4, 2019. Commented out values are from the old setup (just 2016 and 2017).
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016preVFP") 
        {
            norm = 0.04176*nJets + 0.8915;
            expo = (-0.02358*nJets - 0.08940)/1000;
        }
        else if( runYear == "2016postVFP")
        {
            norm = 0.04176*nJets + 0.8915;
            expo = (-0.02358*nJets - 0.08940)/1000;
        }
        else if( runYear == "2017" ) 
        {
            norm = 0.03952*nJets + 0.9171;
            expo = (-0.03275*nJets - 0.03795)/1000;
        }
        else if( runYear == "2018pre" ) 
        {
            norm = 0.03953*nJets + 0.8934;
            expo = (-0.03832*nJets + 0.01108)/1000;
        }
        else if( runYear == "2018post" ) 
        {
            norm = 0.01953*nJets + 1.0150;
            expo = (-0.01409*nJets - 0.1325)/1000;
        }
        return norm*exp( expo*HT ); 
    }

    double htScaleFactorFlat2000(const int nJets, const double HT, const std::string& runYear) const
    {
        return ( HT > 2000.0 ) ? htScaleFactor(nJets,2000.0,runYear) : htScaleFactor(nJets,HT,runYear);
    }

    double htScaleFactorNJet8(const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;      
        if( runYear == "2016preVFP" ) 
        {
            norm = 1.318;
            expo = -0.000349;
        }
        else if( runYear == "2016postVFP" ) 
        {
            norm = 1.318;
            expo = -0.000349;
        }
        else if( runYear == "2017" )
        {
            norm = 1.216;
            expo = -0.0003109;
        }
        else if( runYear == "2018pre" )
        {
            norm = 1.229;
            expo = -0.0002752;
        }
        else if( runYear == "2018post" )
        {
            norm = 1.299;
            expo = -0.0003484;
        }
        return norm*exp( expo*HT );
    }

    double htScaleFactorMG(const int nJets, const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016preVFP" ) 
        {
            norm = -0.002429*nJets + 1.044;
            expo = (0.01214*nJets - 0.1198)/1000;
        }
        else if( runYear == "2016postVFP" ) 
        {
            norm = -0.002429*nJets + 1.044;
            expo = (0.01214*nJets - 0.1198)/1000;
        }
        else if( runYear == "2017" ) 
        {
            norm = 0.1086*nJets + 0.5567;
            expo = (-0.8978*nJets + 0.2198)/1000;
        }
        else if( runYear == "2018pre" ) 
        {
            norm = 0.03349*nJets + 0.9464;
            expo = (-0.02058*nJets - 0.1488)/1000;
        }
        else if( runYear == "2018post" ) 
        {
            norm = 0.02729*nJets + 0.9632;
            expo = (-0.01304*nJets - 0.1637)/1000;
        }
        return norm*exp( expo*HT );
    }

    double htScaleFactorNJet8MG(const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016preVFP" ) 
        {
            norm = 1.073;
            expo = -0.00009169;
        }
        else if( runYear == "2016postVFP" ) 
        {
            norm = 1.073;
            expo = -0.00009169;
        }
        else if( runYear == "2017" )
        {
            norm = 1.236;
            expo = -0.0002949;
        }
        else if( runYear == "2018pre" )
        {
            norm = 1.2;
            expo = -0.0002516;
        }
        else if( runYear == "2018post" )
        {
            norm = 1.273;
            expo = -0.0003242;
        }
        return norm*exp( expo*HT );
    }

    double getMean(const std::string& name)
    {
        if( sfMeanMap_[name] == 0.0 )
        {
            printMeanError_ = true;
            return 1.0;
        }
        else
        {
            return sfMeanMap_[name];
        }
    }

    void scaleFactors(NTupleReader& tr)
    {
        // -----------------------------------------------------------------------------------------------------------------
        // Calculate scale SF and variation
        // Following the example in SusyAnaTools PDFUncertainty.h
        // The scale weights are calculated using the envelope method and we ignore all anti-correlated variations (5 and 7)
        // -----------------------------------------------------------------------------------------------------------------
        const auto& scaleWeights = tr.getVec<float>("ScaleWeights" );
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& runYear      = tr.getVar<std::string>("runYear");
        
        double scaleWeightNominal = 1.0;
        std::vector<float>  myScaleWeights(6, 1.0);
        if( scaleWeights.size() == 9 ) 
        {
            // If there are not exactly 9 scale factors, then the vector  was filled incorrectly
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
            scaleWeightNominal    = 1.0;
            scaleWeightUpperBound = 1.0;
            scaleWeightLowerBound = 1.0;
        }

        double scaleWeightUpperBound_corr = scaleWeightUpperBound;
        double scaleWeightLowerBound_corr = scaleWeightLowerBound;
        if(sfMeanMap_.find(filetag+"_sclUp") != sfMeanMap_.end() && sfMeanMap_.find(filetag+"_sclDown") != sfMeanMap_.end()) 
        {            
            const double meanUp        = getMean(filetag+"_sclUp");
            const double meanDown      = getMean(filetag+"_sclDown");
            scaleWeightUpperBound_corr = (1/meanUp)*scaleWeightUpperBound;
            scaleWeightLowerBound_corr = (1/meanDown)*scaleWeightLowerBound;
        }

        tr.registerDerivedVar("scaleWeightUp"        +myVarSuffix_, scaleWeightUpperBound_corr);
        tr.registerDerivedVar("scaleWeightDown"      +myVarSuffix_, scaleWeightLowerBound_corr);
        tr.registerDerivedVar("scaleWeightUpUncor"   +myVarSuffix_, scaleWeightUpperBound     );
        tr.registerDerivedVar("scaleWeightDownUncor" +myVarSuffix_, scaleWeightLowerBound     );
        tr.registerDerivedVar("scaleWeightNom"       +myVarSuffix_, scaleWeightNominal        );

        // ------------------------------------------------------------------------------
        // Calculate parton shower variation
        // Note: not all samples have these weights stored, give them default value of 1.
        // ------------------------------------------------------------------------------        
        double PSweight_ISRUp   = 1.0, PSweight_ISRDown   = 1.0; 
        double PSweight_FSRUp   = 1.0, PSweight_FSRDown   = 1.0; 
        double PSweight_ISRUp_2 = 1.0, PSweight_ISRDown_2 = 1.0; 
        double PSweight_FSRUp_2 = 1.0, PSweight_FSRDown_2 = 1.0; 
        if(tr.hasVar("PSweights"))
        {
            const auto& PSweights = tr.getVec<float>("PSweights");
            if(PSweights.size() >= 12) // should have size of 14, but just put 12 or more to be able to use the sample with the bug
            {
                // Get nominal one so we can normalize it
                double MEweight = PSweights[0];
                // reduced variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/sqrt(2) and sqrt(2)
                PSweight_ISRUp   = PSweights[2]/MEweight;
                PSweight_FSRUp   = (PSweights[3]/MEweight < 10.0) ? PSweights[3]/MEweight : 1.0;
                PSweight_ISRDown = PSweights[4]/MEweight;
                PSweight_FSRDown = (PSweights[5]/MEweight < 10.0) ? PSweights[5]/MEweight : 1.0;
                // nominal variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/2 and 2
                PSweight_ISRUp_2 = PSweights[6]/MEweight;
                PSweight_FSRUp_2 = (PSweights[7]/MEweight < 10.0) ? PSweights[7]/MEweight : 1.0;
                PSweight_ISRDown_2 = PSweights[8]/MEweight;
                PSweight_FSRDown_2 = (PSweights[9]/MEweight < 10.0) ? PSweights[9]/MEweight : 1.0;
            }
        }
        tr.registerDerivedVar("PSweight_ISRUp"     +myVarSuffix_, PSweight_ISRUp    );
        tr.registerDerivedVar("PSweight_ISRDown"   +myVarSuffix_, PSweight_ISRDown  );
        tr.registerDerivedVar("PSweight_FSRUp"     +myVarSuffix_, PSweight_FSRUp    );
        tr.registerDerivedVar("PSweight_FSRDown"   +myVarSuffix_, PSweight_FSRDown  );
        tr.registerDerivedVar("PSweight_ISRUp_2"   +myVarSuffix_, PSweight_ISRUp_2  );
        tr.registerDerivedVar("PSweight_ISRDown_2" +myVarSuffix_, PSweight_ISRDown_2);
        tr.registerDerivedVar("PSweight_FSRUp_2"   +myVarSuffix_, PSweight_FSRUp_2  );
        tr.registerDerivedVar("PSweight_FSRDown_2" +myVarSuffix_, PSweight_FSRDown_2);

        // ---------------------------------------------------------------------------------------------------
        // Now calculate the PDF scale factor and uncertainty 
        // Based on the 100 different replica values stored in PDFweights using envelope method and the median
        // ---------------------------------------------------------------------------------------------------
        const auto& PDFweights = tr.getVec<float>("PDFweights");
        double central = 1.0, NNPDF_from_median_up = 1.0, NNPDF_from_median_down = 1.0;
        if(PDFweights.size() > 0)
        {
            const double reqCL = 0.68;                        // Choose a confidence level for the uncertainty
            std::vector<float> sortedPDFWeights = PDFweights; // Cannot sort a constant
            std::sort( sortedPDFWeights.begin() + 1, sortedPDFWeights.end() );
        
            const int upper      = std::round( 0.5*(1 + reqCL)*100.0 );
            central              = 0.5*( sortedPDFWeights[50] + sortedPDFWeights[51] ); // Exactly 100 entries
            const double errplus = abs(central - sortedPDFWeights[upper]);

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
            const double meanUp         = getMean(filetag+"_pdf_Up"  );
            const double meanDown       = getMean(filetag+"_pdf_Down");
            NNPDF_from_median_up_corr   = (1/meanUp)*NNPDF_from_median_up;
            NNPDF_from_median_down_corr = (1/meanDown)*NNPDF_from_median_down;
        }

        tr.registerDerivedVar("PDFweightUp"        +myVarSuffix_, NNPDF_from_median_up_corr  );
        tr.registerDerivedVar("PDFweightDown"      +myVarSuffix_, NNPDF_from_median_down_corr);
        tr.registerDerivedVar("PDFweightUpUncor"   +myVarSuffix_, NNPDF_from_median_up       );
        tr.registerDerivedVar("PDFweightDownUncor" +myVarSuffix_, NNPDF_from_median_down     );
        tr.registerDerivedVar("PDFweightNom"       +myVarSuffix_, central                    );
       
        // --------------------------------------------------
        // Now calculate the jet scale factor and uncertainty
        // --------------------------------------------------
        const auto& Jets            = tr.getVec<utility::LorentzVector>(("Jets" +myVarSuffix_));
        const auto& GoodJets_pt45   = tr.getVec<bool>("GoodJets_pt45"           +myVarSuffix_ ); 
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45"          +myVarSuffix_ );
        const auto& HT_trigger_pt45 = tr.getVar<double>("HT_trigger_pt45"       +myVarSuffix_ );


        // find the 6th jet and get its pt
        int njetspt45       = 0;
        double SixthJetPt45 = 0.0;

        for (unsigned int j = 0; j < Jets.size(); j++)
        {
            if (!GoodJets_pt45[j]) continue;
            njetspt45++;

            if (njetspt45 == 6)
            {
                SixthJetPt45 = Jets.at(j).Pt();
                break;
            }
        }
      
        int xbinJetTrig  = 0,   ybinJetTrig   = 0;
        double jetTrigSF = 0.0, jetTrigSF_Err = 0.0; 
        if (NGoodBJets_pt45 == 2)
        {
            xbinJetTrig   = findBin(jetSFHistoTrigName_2bCut_, HT_trigger_pt45, "X", "jet trigger x");
            ybinJetTrig   = findBin(jetSFHistoTrigName_2bCut_, SixthJetPt45,    "Y", "jet trigger y");
            jetTrigSF     = jetSFHistoTrigName_2bCut_->GetBinContent(xbinJetTrig, ybinJetTrig);   
            jetTrigSF_Err = jetSFHistoTrigName_2bCut_->GetBinError(xbinJetTrig, ybinJetTrig); 
        }
       
        else if (NGoodBJets_pt45 == 3)
        {
            xbinJetTrig   = findBin(jetSFHistoTrigName_3bCut_, HT_trigger_pt45, "X", "jet trigger x");
            ybinJetTrig   = findBin(jetSFHistoTrigName_3bCut_, SixthJetPt45,    "Y", "jet trigger y");
            jetTrigSF     = jetSFHistoTrigName_3bCut_->GetBinContent(xbinJetTrig, ybinJetTrig);
            jetTrigSF_Err = jetSFHistoTrigName_3bCut_->GetBinError(xbinJetTrig, ybinJetTrig);
        
        } 

        else if (NGoodBJets_pt45 >= 4)
        {
            xbinJetTrig   = findBin(jetSFHistoTrigName_ge4bCut_, HT_trigger_pt45, "X", "jet trigger x");
            ybinJetTrig   = findBin(jetSFHistoTrigName_ge4bCut_, SixthJetPt45,    "Y", "jet trigger y");
            jetTrigSF     = jetSFHistoTrigName_ge4bCut_->GetBinContent(xbinJetTrig, ybinJetTrig);
            jetTrigSF_Err = jetSFHistoTrigName_ge4bCut_->GetBinError(xbinJetTrig, ybinJetTrig);
        }

        double jetTrigSF_Up   = jetTrigSF + jetTrigSF_Err;
        double jetTrigSF_Down = jetTrigSF - jetTrigSF_Err;

        tr.registerDerivedVar("jetTrigSF"      +myVarSuffix_, jetTrigSF     );
        tr.registerDerivedVar("jetTrigSF_Err"  +myVarSuffix_, jetTrigSF_Err );
        tr.registerDerivedVar("jetTrigSF_Up"   +myVarSuffix_, jetTrigSF_Up  );
        tr.registerDerivedVar("jetTrigSF_Down" +myVarSuffix_, jetTrigSF_Down);

        // -------------------------------------------------------
        // Now calculate the electron scale factor and uncertainty
        // -------------------------------------------------------
        const auto& electrons       = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& goodElectrons   = tr.getVec<bool>("GoodElectrons" +myVarSuffix_);
        double totGoodElectronSF    = 1.0, totGoodElectronSF_Up    = 1.0, totGoodElectronSF_Down    = 1.0;
        double totGoodElectronSFErr = 0.0, totGoodElectronSFPErr2  = 0.0;
        double noTrigGoodElectronSF = 1.0, noTrigGoodElectronSFErr = 0.0, noTrigGoodElectronSFPErr2 = 0.0;
        
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) 
        {
            if( !goodElectrons.at(iel) ) continue;
            
            // Get the scale factor from the rootfile
            const double elpt     = electrons.at(iel).Pt();
            const double eleta    = electrons.at(iel).Eta();
            const int xbinElTight = findBin(eleSFHistoTight_, eleta, "X", "el id x");
            const int ybinElTight = findBin(eleSFHistoTight_, elpt,  "Y", "el id y");                
            const int xbinElIso   = findBin(eleSFHistoIso_,   eleta, "X", "el iso x");
            const int ybinElIso   = findBin(eleSFHistoIso_,   elpt,  "Y", "el iso y");
            const int xbinElReco  = findBin(eleSFHistoReco_,  eleta, "X", "el reco x");
            const int ybinElReco  = findBin(eleSFHistoReco_,  elpt,  "Y", "el reco y");
            const int xbinElTrig  = findBin(eleSFHistoTrig_,  elpt,  "X", "el trigger x");
            const int ybinElTrig  = findBin(eleSFHistoTrig_,  eleta, "Y", "el trigger y");
        
            if( xbinElTight != -1 && ybinElTight != -1 && xbinElIso != -1 && ybinElIso != -1 && xbinElReco != -1 && ybinElReco != -1 ) 
            {
                const double eleTightSF    = eleSFHistoTight_->GetBinContent( xbinElTight, ybinElTight );
                const double eleTightSFErr = eleSFHistoTight_->GetBinError( xbinElTight, ybinElTight );
                const double eleTightPErr  = eleTightSFErr/eleTightSF;
                double eleIsoSF            = eleSFHistoIso_->GetBinContent( xbinElIso, ybinElIso );
                double eleIsoSFErr         = eleSFHistoIso_->GetBinError( xbinElIso, ybinElIso );
                double eleIsoPErr          = eleIsoSFErr/eleIsoSF;
                const double eleRecoSF     = eleSFHistoReco_->GetBinContent( xbinElReco, ybinElReco );
                const double eleRecoSFErr  = eleSFHistoReco_->GetBinError( xbinElReco, ybinElReco );
                const double eleRecoPErr   = eleRecoSFErr/eleRecoSF;
                const double eleTrigSF     = (eleSFHistoTrig_) ? eleSFHistoTrig_->GetBinContent( xbinElTrig, ybinElTrig ) : 1.0;
                const double eleTrigSFErr  = (eleSFHistoTrig_) ? eleSFHistoTrig_->GetBinError( xbinElTrig, ybinElTrig ) : 0.0;
                const double eleTrigPErr   = eleTrigSFErr/eleTrigSF;
                
                if( runYear.find("2016") != std::string::npos ) 
                { 
                    // The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    // If this is the year 2016, we need to add the IP2D histogram scale factors into the Iso scale factor
                    const double eleIP2DSF    = eleSFHistoIP2D_->GetBinContent( xbinElIso, ybinElIso );
                    const double eleIP2DSFErr = eleSFHistoIP2D_->GetBinError( xbinElIso, ybinElIso );
                    const double eleIP2DPErr  = eleIP2DSFErr/eleIP2DSF;

                    eleIsoSF    = eleIsoSF*eleIP2DSF;
                    eleIsoPErr  = utility::addInQuad( eleIsoPErr, eleIP2DPErr );
                    eleIsoSFErr = eleIsoPErr*eleIsoSF;
                }
                
                const double eleNoTrigSF   = eleTightSF*eleIsoSF*eleRecoSF;
                const double eleTotSF      = eleNoTrigSF*eleTrigSF; 
                const double eleNoTrigPErr = utility::addInQuad( eleTightPErr, eleIsoPErr, eleRecoPErr );
                const double eleTotPErr    = utility::addInQuad( eleNoTrigPErr, eleTrigPErr );

                totGoodElectronSF         *= eleTotSF; 
                noTrigGoodElectronSF      *= eleNoTrigSF;
                totGoodElectronSFPErr2    += eleTotPErr*eleTotPErr;
                noTrigGoodElectronSFPErr2 += eleNoTrigPErr*eleNoTrigPErr;
            }            
        }

        totGoodElectronSFErr    = sqrt(totGoodElectronSFPErr2) * totGoodElectronSF;
        noTrigGoodElectronSFErr = sqrt(noTrigGoodElectronSFPErr2) * noTrigGoodElectronSF;
        totGoodElectronSF_Up    = totGoodElectronSF + totGoodElectronSFErr;
        totGoodElectronSF_Down  = totGoodElectronSF - totGoodElectronSFErr;

        tr.registerDerivedVar("totGoodElectronSF"       +myVarSuffix_, totGoodElectronSF      );
        tr.registerDerivedVar("totGoodElectronSFErr"    +myVarSuffix_, totGoodElectronSFErr   );
        tr.registerDerivedVar("totGoodElectronSF_Up"    +myVarSuffix_, totGoodElectronSF_Up   );
        tr.registerDerivedVar("totGoodElectronSF_Down"  +myVarSuffix_, totGoodElectronSF_Down );
        tr.registerDerivedVar("noTrigGoodElectronSF"    +myVarSuffix_, noTrigGoodElectronSF   );
        tr.registerDerivedVar("noTrigGoodElectronSFErr" +myVarSuffix_, noTrigGoodElectronSFErr);

        // -----------------------------------------------
        // Adding code for implementing muon scale factors
        // -----------------------------------------------
        const auto& muons         = tr.getVec<utility::LorentzVector>("Muons" );
        const auto& goodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_  );
        const auto& nonisoMuons   = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);
        double totGoodMuonSF      = 1.0, totGoodMuonSF_Up     = 1.0, totGoodMuonSF_Down    = 1.0;
        double totGoodMuonSFErr   = 0.0, totGoodMuonSFPErr2   = 0.0;
        double noTrigGoodMuonSF   = 1.0, noTrigGoodMuonSFErr  = 0.0, noTrigGoodMuonSFPErr2 = 0.0;
        double totNonIsoMuonSF    = 1.0, totNonIsoMuonSF_Up   = 1.0, totNonIsoMuonSF_Down  = 1.0;
        double totNonIsoMuonSFErr = 0.0, totNonIsoMuonSFPErr2 = 0.0;        

        for( unsigned int imu = 0; imu < muons.size(); imu++ ) 
        {            
            // Get the scale factor from the rootfile
            const double mupt = muons.at(imu).Pt();
            const double mueta = muons.at(imu).Eta();

            const int ybinMuMedium     = findBin(muSFHistoMedium_, mupt,       "Y", "mu id y"         );
            const int xbinMuMedium     = findBin(muSFHistoMedium_, abs(mueta), "X", "mu id x"         );
            const int ybinMuIso        = findBin(muSFHistoIso_,    mupt,       "Y", "mu iso y"        );
            const int xbinMuIso        = findBin(muSFHistoIso_,    abs(mueta), "X", "mu iso x"        );            
            const int ybinMuTrig       = findBin(muSFHistoTrig_,   mueta,      "Y", "mu trigger y"    );
            const int xbinMuTrig       = findBin(muSFHistoTrig_,   mupt,       "X", "mu trigger x"    );
            const int ybinNonIsoMuTrig = findBin(nimuSFHistoTrig_, mueta,      "Y", "mu iso trigger y");
            const int xbinNonIsoMuTrig = findBin(nimuSFHistoTrig_, mupt,       "X", "mu iso trigger x");
            if( xbinMuMedium != -1 && ybinMuMedium != -1 && xbinMuIso != -1 && ybinMuIso != -1 ) 
            {
                // The SUSLepton Twiki claims that the errors in the histogrm are purely statistical and can be ignored and recommends a 3% error for each leg (ID+IP+ISO)
                const double muMediumSF         = muSFHistoMedium_->GetBinContent( xbinMuMedium, ybinMuMedium  );
                const double muIsoSF            = muSFHistoIso_->GetBinContent( xbinMuIso, ybinMuIso           );
                const double muTrigSF           = (muSFHistoTrig_) ? muSFHistoTrig_->GetBinContent( xbinMuTrig, ybinMuTrig ) : 1.0;
                const double muNonIsoTrigSF     = (nimuSFHistoTrig_) ? nimuSFHistoTrig_->GetBinContent( xbinNonIsoMuTrig, ybinNonIsoMuTrig ) : 1.0;
                const double muTrigSFErr        = (muSFHistoTrig_) ? muSFHistoTrig_->GetBinError( xbinMuTrig, ybinMuTrig ) : 0.0;
                const double muNonIsoTrigSFErr  = (nimuSFHistoTrig_) ? nimuSFHistoTrig_->GetBinError( xbinNonIsoMuTrig, ybinNonIsoMuTrig ) : 0.0;
                const double muTrigSFPErr       = muTrigSFErr/muTrigSF;
                const double muNonIsoTrigSFPErr = muNonIsoTrigSFErr/muNonIsoTrigSF;
                double muNoTrigSF               = muMediumSF*muIsoSF; 
                double muTotSF                  = muNoTrigSF*muTrigSF;
                double muNonIsoTotSF            = muNoTrigSF*muNonIsoTrigSF;
                const double muNoTrigSFPErr2    = 0.03*0.03;
                const double muTotSFPErr2       = muNoTrigSFPErr2 + muTrigSFPErr*muTrigSFPErr;
                const double muNonIsoTotSFPErr2 = muNoTrigSFPErr2 + muNonIsoTrigSFPErr*muNonIsoTrigSFPErr;
                
                
                // The RECO scale factor for muons have been shown to be close to unity and are not required by the Muon POG for UL (see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonTnPOverview). Confirmation from the SUSY Muon contact is still needed

                /*if( runYear.find("2016") != std::string::npos ) 
                {
                    //For the general track reconstruction they claim that the errors for the systematic still need to be finalized - does not seem to have been finalized as of Dec 2018
                    //This reconstruction value only exists for 2016 - SUS SF people say the 3% will include the reco scale factor uncertainty for now
                    const double muRecoSF    = muSFHistoReco_->Eval( mueta );
                    muTotSF                 *= muRecoSF;
                    muNoTrigSF              *= muRecoSF;
                    muNonIsoTotSF           *= muRecoSF;
                }
                */
                
                if( goodMuons.at(imu) )
                {                
                    totGoodMuonSF         *= muTotSF;
                    noTrigGoodMuonSF      *= muNoTrigSF;
                    totGoodMuonSFPErr2    += muTotSFPErr2;
                    noTrigGoodMuonSFPErr2 += muNoTrigSFPErr2;
                }
                if( nonisoMuons.at(imu) )
                {
                    totNonIsoMuonSF      *= muNonIsoTotSF;
                    totNonIsoMuonSFPErr2 += muNonIsoTotSFPErr2;
                }
            }
        }

        totGoodMuonSFErr     = sqrt(totGoodMuonSFPErr2)*totGoodMuonSF;
        totGoodMuonSF_Up     = totGoodMuonSF + totGoodMuonSFErr;
        totGoodMuonSF_Down   = totGoodMuonSF - totGoodMuonSFErr;
        noTrigGoodMuonSFErr  = sqrt(noTrigGoodMuonSFPErr2)*noTrigGoodMuonSF;
        totNonIsoMuonSFErr   = sqrt(totNonIsoMuonSFPErr2)*totNonIsoMuonSF;
        totNonIsoMuonSF_Up   = totNonIsoMuonSF + totNonIsoMuonSFErr;
        totNonIsoMuonSF_Down = totNonIsoMuonSF - totNonIsoMuonSFErr;

        tr.registerDerivedVar("totGoodMuonSF"        +myVarSuffix_, totGoodMuonSF       );
        tr.registerDerivedVar("totGoodMuonSFErr"     +myVarSuffix_, totGoodMuonSFErr    );
        tr.registerDerivedVar("totGoodMuonSF_Up"     +myVarSuffix_, totGoodMuonSF_Up    );
        tr.registerDerivedVar("totGoodMuonSF_Down"   +myVarSuffix_, totGoodMuonSF_Down  );
        tr.registerDerivedVar("noTrigGoodMuonSF"     +myVarSuffix_, noTrigGoodMuonSF    );
        tr.registerDerivedVar("noTrigGoodMuonSFErr"  +myVarSuffix_, noTrigGoodMuonSFErr );
        tr.registerDerivedVar("totNonIsoMuonSF"      +myVarSuffix_, totNonIsoMuonSF     );
        tr.registerDerivedVar("totNonIsoMuonSFErr"   +myVarSuffix_, totNonIsoMuonSFErr  );
        tr.registerDerivedVar("totNonIsoMuonSF_Up"   +myVarSuffix_, totNonIsoMuonSF_Up  );
        tr.registerDerivedVar("totNonIsoMuonSF_Down" +myVarSuffix_, totNonIsoMuonSF_Down);

        // -------------------------------------------------------------------------------
        // Adding a scale factor that corrects the disagreement between data and MC for Ht
        // -------------------------------------------------------------------------------
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& isSignal        = tr.getVar<bool>("isSignal");

        // Get the uncorrected ht scale factors
        const double htDerivedweightUncor         = htScaleFactor(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double htDerivedweightFlat2000Uncor = htScaleFactorFlat2000(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double htDerivedweightNJet7Uncor    = htScaleFactor(7, HT_trigger_pt30, runYear);
        const double htDerivedweightMGUncor       = htScaleFactorMG(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double fit2NJetBin8                 = htScaleFactorNJet8(HT_trigger_pt30, runYear);
        const double fit2NJetBin567               = htScaleFactor(8, HT_trigger_pt30, runYear);
        const double fit2NJetBin8MG               = htScaleFactorNJet8MG(HT_trigger_pt30, runYear);
        const double fit2NJetBin567MG             = htScaleFactor(8, HT_trigger_pt30, runYear);
                
        double htDerivedweight = 1.0, htDerivedweightFlat2000 = 1.0, htDerivedweightNJet7 = 1.0, htDerivedweightMG = 1.0;
        double htScaleUpUncor = 1.0, htScaleDownUncor = 1.0, htScaleUpMGUncor = 1.0, htScaleDownMGUncor = 1.0;
        double htScaleUp = 1.0, htScaleDown = 1.0, htScaleUpMG = 1.0, htScaleDownMG = 1.0;
        if( sfMeanMap_.find(filetag+"_ht") != sfMeanMap_.end() && !isSignal ) 
        {
            htDerivedweight         = (1.0/getMean(filetag+"_ht"         ))*htDerivedweightUncor;
            htDerivedweightFlat2000 = (1.0/getMean(filetag+"_ht_flat2000"))*htDerivedweightFlat2000Uncor;
            htDerivedweightNJet7    = (1.0/getMean(filetag+"_ht_njet7"   ))*htDerivedweightNJet7Uncor;
            htDerivedweightMG       = (1.0/getMean(filetag+"_ht_MG"      ))*htDerivedweightMGUncor;
            htScaleUpUncor          = htDerivedweight*(fit2NJetBin8/fit2NJetBin567);
            htScaleDownUncor        = htDerivedweight*(fit2NJetBin567/fit2NJetBin8);
            htScaleUpMGUncor        = htDerivedweightMG*(fit2NJetBin8MG/fit2NJetBin567MG);
            htScaleDownMGUncor      = htDerivedweightMG*(fit2NJetBin567MG/fit2NJetBin8MG);
            htScaleUp               = (1.0/getMean(filetag+"_htUp"     ))*htScaleUpUncor;
            htScaleDown             = (1.0/getMean(filetag+"_htDown"   ))*htScaleDownUncor;
            htScaleUpMG             = (1.0/getMean(filetag+"_ht_MGUp"  ))*htScaleUpMGUncor;
            htScaleDownMG           = (1.0/getMean(filetag+"_ht_MGDown"))*htScaleDownMGUncor;
        }

        tr.registerDerivedVar("htDerivedweight"              +myVarSuffix_, htDerivedweight             );
        tr.registerDerivedVar("htDerivedweightFlat2000"      +myVarSuffix_, htDerivedweightFlat2000     );
        tr.registerDerivedVar("htDerivedweightNJet7"         +myVarSuffix_, htDerivedweightNJet7        );
        tr.registerDerivedVar("htDerivedweightMG"            +myVarSuffix_, htDerivedweightMG           );
        tr.registerDerivedVar("htDerivedweightUncor"         +myVarSuffix_, htDerivedweightUncor        );
        tr.registerDerivedVar("htDerivedweightFlat2000Uncor" +myVarSuffix_, htDerivedweightFlat2000Uncor);
        tr.registerDerivedVar("htDerivedweightNJet7Uncor"    +myVarSuffix_, htDerivedweightNJet7Uncor   );
        tr.registerDerivedVar("htDerivedweightMGUncor"       +myVarSuffix_, htDerivedweightMGUncor      );
        tr.registerDerivedVar("htScaleUpUncor"               +myVarSuffix_, htScaleUpUncor              );
        tr.registerDerivedVar("htScaleDownUncor"             +myVarSuffix_, htScaleDownUncor            );
        tr.registerDerivedVar("htScaleUpMGUncor"             +myVarSuffix_, htScaleUpMGUncor            );
        tr.registerDerivedVar("htScaleDownMGUncor"           +myVarSuffix_, htScaleDownMGUncor          );
        tr.registerDerivedVar("htScaleUp"                    +myVarSuffix_, htScaleUp                   );
        tr.registerDerivedVar("htScaleDown"                  +myVarSuffix_, htScaleDown                 );
        tr.registerDerivedVar("htScaleUpMG"                  +myVarSuffix_, htScaleUpMG                 );
        tr.registerDerivedVar("htScaleDownMG"                +myVarSuffix_, htScaleDownMG               );

        //---------------------------------------------------------------------------------------------------------
        // Adding a scale factor for pileup 
        // For 2016 and 2018: Grab the individual pileup weight from the histogram found in PileupHistograms_*.root
        // For 2017: Using the puWeight stored in the nTuples
        // --------------------------------------------------------------------------------------------------------        
        const auto& puWeightUnCorr  = tr.getVar<float>("puWeight" );
        const auto& puSysUpUnCorr   = tr.getVar<float>("puSysUp"  );
        const auto& puSysDownUnCorr = tr.getVar<float>("puSysDown");

        tr.registerDerivedVar( "puWeightUnCorr" +myVarSuffix_, puWeightUnCorr );
        tr.registerDerivedVar( "puSysUpUnCorr"  +myVarSuffix_, puSysUpUnCorr  );
        tr.registerDerivedVar( "puSysDownUnCorr"+myVarSuffix_, puSysDownUnCorr);

        // Adding correction to the pileup weight
        double puWeightCorr  = puWeightUnCorr;
        double puSysUpCorr   = puSysUpUnCorr;
        double puSysDownCorr = puSysDownUnCorr;
        if( sfMeanMap_.find(filetag+"_pu") != sfMeanMap_.end() && !isSignal )
        {
            const double mean     = getMean(filetag+"_pu"     );
            const double meanUp   = getMean(filetag+"_pu_Up"  );
            const double meanDown = getMean(filetag+"_pu_Down");
            puWeightCorr = (1.0/mean)*puWeightUnCorr;
            puSysUpCorr  = (1.0/meanUp)*puSysUpUnCorr;
            puSysDownCorr= (1.0/meanDown)*puSysDownUnCorr;
        }
        tr.registerDerivedVar("puWeightCorr"  +myVarSuffix_, puWeightCorr );
        tr.registerDerivedVar("puSysUpCorr"   +myVarSuffix_, puSysUpCorr  );
        tr.registerDerivedVar("puSysDownCorr" +myVarSuffix_, puSysDownCorr);
        
        // -------------------------------------------------------------
        // Adding top pt reweighting for ttbar MC (Powheg)
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        // 13 TeV all combined
        // -------------------------------------------------------------
        double topPtScaleFactor = 1.0;
        auto* topPtVec          = new std::vector<float>();
        if( filetag.find("TTTo") != std::string::npos )
        {
            const double a=0.103, b=-0.0118, c=-0.000134, d=0.973;
            auto SF = [&](const double pt){return a * exp(b * pt) + c * pt + d;};
            
            const auto& GenParticles        = tr.getVec<utility::LorentzVector>("GenParticles");
            const auto& GenParticles_PdgId  = tr.getVec<int>("GenParticles_PdgId"             );
            const auto& GenParticles_Status = tr.getVec<int>("GenParticles_Status"            );

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++)
            {   
                if( abs(GenParticles_PdgId[gpi]) == 6 && GenParticles_Status[gpi] == 62 )
                {
                    topPtScaleFactor *= SF( GenParticles[gpi].Pt() );
                    topPtVec->push_back( GenParticles[gpi].Pt() );
                }
            }
            topPtScaleFactor = sqrt(topPtScaleFactor);
        }
        
        tr.registerDerivedVar("topPtScaleFactor" +myVarSuffix_, topPtScaleFactor);
        tr.registerDerivedVec("topPtVec"         +myVarSuffix_, topPtVec        );
       
        // -----------------------------------------------------------------------------------------------------------------
        // Adding reweighting recipe to emulate Level 1 ECAL prefiring
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
        // All prefiring recipes are now included in UL ntuples. Values are grabbed directly and registered with tree reader
        // -----------------------------------------------------------------------------------------------------------------
        double prefiringScaleFactor = 1.0, prefiringScaleFactorUp = 1.0, prefiringScaleFactorDown = 1.0;
        prefiringScaleFactor        = tr.getVar<float>("NonPrefiringProb"    );
        prefiringScaleFactorUp      = tr.getVar<float>("NonPrefiringProbUp"  );
        prefiringScaleFactorDown    = tr.getVar<float>("NonPrefiringProbDown");
        tr.registerDerivedVar("prefiringScaleFactor"     +myVarSuffix_, prefiringScaleFactor    );                    
        tr.registerDerivedVar("prefiringScaleFactorUp"   +myVarSuffix_, prefiringScaleFactorUp  );                    
        tr.registerDerivedVar("prefiringScaleFactorDown" +myVarSuffix_, prefiringScaleFactorDown);                    

        const auto& Weight            = tr.getVar<float>("Weight");
        const auto& FinalLumi         = tr.getVar<double>("FinalLumi");

        const auto& analyzer = tr.getVar<std::string>("analyzer");

        // For the CalculateSFMean analyzer, we do not need to run BTagCorrector
        // and thus, we do not need to try and get the bTagSF here, so just default 1.0
        double bTagWeight = 1.0;
        if ( analyzer != "CalculateSFMean" )
        {
            bTagWeight = tr.getVar<double>("bTagSF_EventWeightSimple_Central" +myVarSuffix_);
        }

        const auto& NGoodMuons        = tr.getVar<int>("NGoodMuons"                          +myVarSuffix_);
        const auto& NGoodLeptons      = tr.getVar<int>("NGoodLeptons"                        +myVarSuffix_);
        const auto& NNonIsoMuons      = tr.getVar<int>("NNonIsoMuons"                        +myVarSuffix_);
        double totalEventWeight       = -1.0;

        double commonWeight = Weight * FinalLumi * bTagWeight * prefiringScaleFactor * puWeightCorr;
        tr.registerDerivedVar("CommonWeight" + myVarSuffix_, commonWeight);

        // 0-Lepton
        // Eventually, jetTrigSF will become totJetSF to incorporate the top tagger scale factor
        if ( NGoodLeptons == 0 and NNonIsoMuons == 0 ) 
        {
            totalEventWeight = commonWeight * jetTrigSF;
        }
        // 1-Lepton
        // Just choose either the muon or the electron weight
        else if ( NGoodLeptons == 1 )
        { 
            double totLepWeight = totGoodElectronSF;
            if ( NGoodMuons == 1 ) 
                totLepWeight = totGoodMuonSF;

            totalEventWeight   = commonWeight * totLepWeight;
        }
        // 2-Lepton
        // TODO: Figure out how to incorporate lep scale factor for more than one lepton
        else if ( NGoodLeptons == 2 )
        {
            totalEventWeight   = commonWeight;// * totLepWeight;
        }
        // QCDCR
        // With one non-isolated muon, use the totNonIsoMuonSF inside the total weight
        else if ( NGoodLeptons == 0 and NNonIsoMuons == 1 ) 
        {
            totalEventWeight    = commonWeight * totNonIsoMuonSF;
        }

        tr.registerDerivedVar("TotalWeight"           + myVarSuffix_, totalEventWeight       );

        if(printGetBinError_) firstPrint_ = false;
    }

public:
    ScaleFactors( const std::string& runYear, const std::string& leptonFileName, const std::string& hadronicFileName, const std::string& meanFileName, const std::string& myVarSuffix = "" )
        : myVarSuffix_(myVarSuffix), firstPrint_(true), printMeanError_(false), printGetBinError_(false)
    {
        std::cout<<"Setting up ScaleFactors"<<std::endl;
        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?

        // Getting Leptonic and Hadronic scale factor histograms
        TFile leptonic_SFRootFile( leptonFileName.c_str()    );
        TFile hadronic_SFRootFile( hadronicFileName.c_str() );
        TString eleSFHistoTightName,      eleSFHistoIsoName,        eleSFHistoRecoName,         eleSFHistoTrigName;
        TString muSFHistoMediumName,      muSFHistoIsoName,         muSFHistoRecoName,          muSFHistoTrigName, nimuSFHistoTrigName;
        TString jetSFHistoTrigName_2bCut, jetSFHistoTrigName_3bCut, jetSFHistoTrigName_ge4bCut;

        // Electron Iso root file not currently seen in recommendations. Working to determine if this is still necessary
        if( runYear == "2016preVFP")
        {
            eleSFHistoTightName        = "EGamma_SF2D_2016preVFP_UL_ID";
            //eleSFHistoIsoName         = "Run2016preVFP_Mini";
            eleSFHistoRecoName         = "EGamma_SF2D_2016preVFP_UL_RECO";
            eleSFHistoTrigName         = "TrigEff_2016preVFP_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName        = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_2016preVFP_UL_ID";
            muSFHistoIsoName           = "NUM_TightRelIso_DEN_MediumID_abseta_pt_2016preVFP_UL_ISO";
            muSFHistoTrigName          = "TrigEff_2016preVFP_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            //nimuSFHistoTrigName       = "TrigEff_2016preVFP_num_nimu_pt40_trig_4jCut";
            nimuSFHistoTrigName        = ""; //just for calculating non iso muon scale factors
            jetSFHistoTrigName_2bCut   = "h_2016preVFP_CombHadIsoMu_trig_2bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_3bCut   = "h_2016preVFP_CombHadIsoMu_trig_3bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_ge4bCut = "h_2016preVFP_CombHadIsoMu_trig_ge4bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT"; 
        }
        else if( runYear == "2016postVFP")
        {
            eleSFHistoTightName        = "EGamma_SF2D_2016postVFP_UL_ID";
            //eleSFHistoIsoName          = "Run2016postVFP_Mini";
            eleSFHistoRecoName         = "EGamma_SF2D_2016postVFP_UL_RECO";
            eleSFHistoTrigName         = "TrigEff_2016postVFP_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName        = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_2016postVFP_UL_ID";
            muSFHistoIsoName           = "NUM_TightRelIso_DEN_MediumID_abseta_pt_2016postVFP_UL_ISO";
            muSFHistoTrigName          = "TrigEff_2016postVFP_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            //nimuSFHistoTrigName        = "TrigEff_2016postVFP_num_nimu_pt40_trig_4jCut";
            nimuSFHistoTrigName        = ""; //just for calculating non iso muon scale factors
            jetSFHistoTrigName_2bCut   = "h_2016postVFP_CombHadIsoMu_trig_2bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_3bCut   = "h_2016postVFP_CombHadIsoMu_trig_3bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_ge4bCut = "h_2016postVFP_CombHadIsoMu_trig_ge4bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
        }
        else if( runYear == "2017")
        {
            eleSFHistoTightName        = "EGamma_SF2D_2017_UL_ID";
            //eleSFHistoIsoName          = "Run2017_MVAVLooseTightIP2DMini";
            eleSFHistoRecoName         = "EGamma_SF2D_2017_UL_RECO";
            eleSFHistoTrigName         = "TrigEff_2017_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName        = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_2017_UL_ID";
            muSFHistoIsoName           = "NUM_TightRelIso_DEN_MediumID_abseta_pt_2017_UL_ISO";
            muSFHistoTrigName          = "TrigEff_2017_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            //nimuSFHistoTrigName        = "TrigEff_2017_num_nimu_pt40_trig_4jCut";
            nimuSFHistoTrigName        = ""; //just for calculating non iso muon scale factors
            jetSFHistoTrigName_2bCut   = "h_2017_CombHadIsoMu_trig_2bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_3bCut   = "h_2017_CombHadIsoMu_trig_3bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_ge4bCut = "h_2017_CombHadIsoMu_trig_ge4bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
        }
        else if( runYear == "2018")
        {
            eleSFHistoTightName        = "EGamma_SF2D_2018_UL_ID";
            //eleSFHistoIsoName          = "Run2018_Mini";
            eleSFHistoRecoName         = "EGamma_SF2D_2018_UL_RECO";
            eleSFHistoTrigName         = "TrigEff_2018_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName        = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_2018_UL_ID";
            muSFHistoIsoName           = "NUM_TightRelIso_DEN_MediumID_abseta_pt_2018_UL_ISO";
            muSFHistoTrigName          = "TrigEff_2018_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            nimuSFHistoTrigName        = "";
            jetSFHistoTrigName_2bCut   = "h_2018_CombHadIsoMu_trig_2bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_3bCut   = "h_2018_CombHadIsoMu_trig_3bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
            jetSFHistoTrigName_ge4bCut = "h_2018_CombHadIsoMu_trig_ge4bjetCut_pt45_HTvs6thJetPt_SingleMuon_TT";
        }

        getHisto(leptonic_SFRootFile, eleSFHistoTight_,            eleSFHistoTightName       );
        //getHisto(leptonic_SFRootFile, eleSFHistoIso_,              eleSFHistoIsoName         );
        getHisto(leptonic_SFRootFile, eleSFHistoReco_,             eleSFHistoRecoName        );
        getHisto(leptonic_SFRootFile, eleSFHistoTrig_,             eleSFHistoTrigName        );
        getHisto(leptonic_SFRootFile, muSFHistoMedium_,            muSFHistoMediumName       );
        getHisto(leptonic_SFRootFile, muSFHistoIso_,               muSFHistoIsoName          );
        getHisto(leptonic_SFRootFile, muSFHistoTrig_,              muSFHistoTrigName         );
        getHisto(leptonic_SFRootFile, nimuSFHistoTrig_,            nimuSFHistoTrigName       );
        getHisto(hadronic_SFRootFile, jetSFHistoTrigName_2bCut_,   jetSFHistoTrigName_2bCut  );
        getHisto(hadronic_SFRootFile, jetSFHistoTrigName_3bCut_,   jetSFHistoTrigName_3bCut  );
        getHisto(hadronic_SFRootFile, jetSFHistoTrigName_ge4bCut_, jetSFHistoTrigName_ge4bCut);        

        leptonic_SFRootFile.Close();
        hadronic_SFRootFile.Close();        

        // Getting mean of some scale factors to keep total number of events the same after apply the scale factor
        TFile SFMeanRootFile( meanFileName.c_str() );
        TIter next(SFMeanRootFile.GetListOfKeys());
        TKey* key;
        while((key = static_cast<TKey*>(next())))
        {
            std::shared_ptr<TH1> h( static_cast<TH1*>(key->ReadObj()) );
            std::string name( h->GetTitle() );
            sfMeanMap_.insert(std::pair<std::string, double>(name, h->GetMean()));
        }
        SFMeanRootFile.Close();
    }

    ~ScaleFactors() 
    {
        if(printMeanError_)   std::cerr<<utility::color("Error: Scale Factor mean is 0.0 setting it to 1.0", "red")<<std::endl;

        std::string message = "";
        for(const auto& m : binError_message_) message += (message=="") ? m : ", "+m;
        if(printGetBinError_) std::cerr<<utility::color("Warning: There was an error in extracting the bin index for the following scale factor(s): \""+message+"\"", "red")<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + myVarSuffix_);
        const auto& fastMode       = tr.getVar<bool>("fastMode");

        if (!lostCauseEvent or !fastMode)
            scaleFactors(tr);
    }
};

#endif
