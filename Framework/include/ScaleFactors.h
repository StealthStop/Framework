#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
#include "TFile.h"
#include "TKey.h"
#include "Framework/Framework/include/Utility.h"

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"

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
    std::shared_ptr<TH2F> jetSFHistoTrigName_1bCut_;
    std::shared_ptr<TH2F> jetSFHistoTrigName_ge2bCut_;
    std::shared_ptr<TH1F> topTagSFHisto_Res_;
    std::shared_ptr<TH1F> topTagSFHisto_Mrg_;
    std::shared_ptr<TH1F> topMistagSFHisto_Res_;
    std::shared_ptr<TH1F> topMistagSFHisto_Mrg_;
    std::shared_ptr<TH2F> topTagEffHisto_Mrg_den_;
    std::shared_ptr<TH2F> topTagEffHisto_Res_den_;
    std::shared_ptr<TH2F> topTagMisHisto_Mrg_den_;
    std::shared_ptr<TH2F> topTagMisHisto_Res_den_;
    std::shared_ptr<TH2F> topTagEffHisto_Mrg_num_;
    std::shared_ptr<TH2F> topTagEffHisto_Res_num_;
    std::shared_ptr<TH2F> topTagMisHisto_Mrg_num_;
    std::shared_ptr<TH2F> topTagMisHisto_Res_num_;

    template<typename T> std::shared_ptr<T>& getHisto(TFile& f, std::shared_ptr<T>& h, const TString& name)
    {
        const auto& ptr = static_cast<T*>(f.Get(name));
        if(name != "" and ptr)
        {
            h.reset( ptr );
        }
        else
        {
            std::cerr<<utility::color("Warning: A needed scale factor histogram, \"" + std::string(name) + "\", is set to nullptr, therefore using 1.0 as the default", "red")<<std::endl;
            h = nullptr;
        }

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
              auto  filetagTemp  = tr.getVar<std::string>("filetag");
        const auto& runYear      = tr.getVar<std::string>("runYear");

        // If we run on a skim data set, there should be a "_skim" suffix
        // at the end of the name, so scrape it off before trying
        // to get anything that requires the filetag
        std::string scrapeOff = "_skim";
        const size_t position = filetagTemp.find(scrapeOff);
        if(position != std::string::npos)
        {
            filetagTemp.erase(position, scrapeOff.length());
        }
        std::string filetag = filetagTemp;

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
        if (NGoodBJets_pt45 == 1)
        {
            xbinJetTrig   = findBin(jetSFHistoTrigName_1bCut_, HT_trigger_pt45, "X", "jet trigger x");
            ybinJetTrig   = findBin(jetSFHistoTrigName_1bCut_, SixthJetPt45,    "Y", "jet trigger y");

            if (xbinJetTrig != -1 and ybinJetTrig != -1)
            {
                jetTrigSF     = jetSFHistoTrigName_1bCut_->GetBinContent(xbinJetTrig, ybinJetTrig       );   
                jetTrigSF_Err = jetSFHistoTrigName_1bCut_->GetBinError(xbinJetTrig, ybinJetTrig         ); 
            }
        }

        else if (NGoodBJets_pt45 >= 2)
        {
            xbinJetTrig   = findBin(jetSFHistoTrigName_ge2bCut_, HT_trigger_pt45, "X", "jet trigger x");
            ybinJetTrig   = findBin(jetSFHistoTrigName_ge2bCut_, SixthJetPt45,    "Y", "jet trigger y");

            if (xbinJetTrig != -1 and ybinJetTrig != -1)
            {
                jetTrigSF     = jetSFHistoTrigName_ge2bCut_->GetBinContent(xbinJetTrig, ybinJetTrig       );
                jetTrigSF_Err = jetSFHistoTrigName_ge2bCut_->GetBinError(xbinJetTrig, ybinJetTrig         );
            }
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
        // Adding a top tagging scale factor in the spirit of b tag sf via "method 1a"
        // -------------------------------------------------------------------------------
        double mcTag     = 1.0, mcNoTag     = 1.0, dataTag     = 1.0, dataNoTag     = 1.0;
        double mcTagUp   = 1.0, mcNoTagUp   = 1.0, dataTagUp   = 1.0, dataNoTagUp   = 1.0;
        double mcTagDown = 1.0, mcNoTagDown = 1.0, dataTagDown = 1.0, dataNoTagDown = 1.0;

        // Distinguish if the best top candidate is resolved or merged
        const auto& resolvedWP = tr.getVar<double>("resolvedTop_WP");
        const auto& mergedWP   = tr.getVar<double>("mergedTop_WP");

        //loop over jets
        const auto* topTagRes = tr.getVar<TopTaggerResults*>("ttr");
        const auto& tops      = topTagRes->getTops();
        for(const auto& top : tops)
        {
            const auto* genTop = top->getBestGenTopMatch();

            bool isResolved = top->getType() == TopObject::RESOLVED_TOP;
            bool isMerged   = top->getType() == TopObject::MERGED_TOP;

            double num = 1.0, numUnc = 0.0, den = 1.0, denUnc = 0.0, sf = 1.0, sfUnc = 0.0;
            int xBinTopNum = -1, yBinTopNum = -1, xBinTopDen = -1, yBinTopDen = -1, binTopSF = -1;
            // Efficiency when dealing with an actual top
            if (genTop)
            {
                if (isResolved)
                {
                    xBinTopNum = findBin(topTagEffHisto_Res_num_, top->P().Pt(),  "X", "top tag eff num x");
                    yBinTopNum = findBin(topTagEffHisto_Res_num_, top->P().Eta(), "Y", "top tag eff num y");
                    xBinTopDen = findBin(topTagEffHisto_Res_den_, top->P().Pt(),  "X", "top tag eff den x");
                    yBinTopDen = findBin(topTagEffHisto_Res_den_, top->P().Eta(), "Y", "top tag eff den y");
                
                    if (xBinTopNum != -1 and yBinTopNum != -1 and xBinTopDen != -1 and yBinTopDen != -1)
                    {
                        num    = topTagEffHisto_Res_num_->GetBinContent(xBinTopNum, yBinTopNum);
                        numUnc = topTagEffHisto_Res_num_->GetBinError(xBinTopNum,   yBinTopNum);
                        den    = topTagEffHisto_Res_den_->GetBinContent(xBinTopDen, yBinTopDen);
                        denUnc = topTagEffHisto_Res_den_->GetBinError(xBinTopDen,   yBinTopDen);
                    }

                    binTopSF = findBin(topTagSFHisto_Res_, top->P().Pt(), "X", "top tag sf x");
                    if (binTopSF != -1)
                    {
                        sf    = topTagSFHisto_Res_->GetBinContent(binTopSF);
                        sfUnc = topTagSFHisto_Res_->GetBinError(binTopSF);
                    }
                }
                else if (isMerged)
                {
                    xBinTopNum = findBin(topTagEffHisto_Mrg_num_, top->P().Pt(),  "X", "top tag eff num x");
                    yBinTopNum = findBin(topTagEffHisto_Mrg_num_, top->P().Eta(), "Y", "top tag eff num y");
                    xBinTopDen = findBin(topTagEffHisto_Mrg_den_, top->P().Pt(),  "X", "top tag eff den x");
                    yBinTopDen = findBin(topTagEffHisto_Mrg_den_, top->P().Eta(), "Y", "top tag eff den y");

                    if (xBinTopNum != -1 and yBinTopNum != -1 and xBinTopDen != -1 and yBinTopDen != -1)
                    {
                        num    = topTagEffHisto_Mrg_num_->GetBinContent(xBinTopNum, yBinTopNum);
                        numUnc = topTagEffHisto_Mrg_num_->GetBinError(xBinTopNum,   yBinTopNum);
                        den    = topTagEffHisto_Mrg_den_->GetBinContent(xBinTopDen, yBinTopDen);
                        denUnc = topTagEffHisto_Mrg_den_->GetBinError(xBinTopDen,   yBinTopDen);
                    }

                    binTopSF = findBin(topTagSFHisto_Mrg_, top->P().Pt(), "X", "top tag sf x");
                    if (binTopSF != -1)
                    {
                        sf    = topTagSFHisto_Mrg_->GetBinContent(binTopSF);
                        sfUnc = topTagSFHisto_Mrg_->GetBinError(binTopSF);
                    }
                }
            }
            // Mistag when dealing with a fake top i.e. no GEN top present
            else
            {
                if (isResolved)
                {
                    xBinTopNum = findBin(topTagMisHisto_Res_num_, top->P().Pt(),  "X", "top tag eff num x");
                    yBinTopNum = findBin(topTagMisHisto_Res_num_, top->P().Eta(), "Y", "top tag eff num y");
                    xBinTopDen = findBin(topTagMisHisto_Res_den_, top->P().Pt(),  "X", "top tag eff den x");
                    yBinTopDen = findBin(topTagMisHisto_Res_den_, top->P().Eta(), "Y", "top tag eff den y");

                    if (xBinTopNum != -1 and yBinTopNum != -1 and xBinTopDen != -1 and yBinTopDen != -1)
                    {
                        num    = topTagMisHisto_Res_num_->GetBinContent(xBinTopNum, yBinTopNum);
                        numUnc = topTagMisHisto_Res_num_->GetBinError(xBinTopNum,   yBinTopNum);
                        den    = topTagMisHisto_Res_den_->GetBinContent(xBinTopDen, yBinTopDen);
                        denUnc = topTagMisHisto_Res_den_->GetBinError(xBinTopDen,   yBinTopDen);
                    }

                    binTopSF = findBin(topMistagSFHisto_Res_, top->P().Pt(), "X", "top tag sf x");
                    if (binTopSF != -1)
                    {
                        sf    = topMistagSFHisto_Res_->GetBinContent(binTopSF);
                        sfUnc = topMistagSFHisto_Res_->GetBinError(binTopSF);
                    }
                }
                else if (isMerged)
                {
                    xBinTopNum = findBin(topTagMisHisto_Mrg_num_, top->P().Pt(),  "X", "top tag eff num x");
                    yBinTopNum = findBin(topTagMisHisto_Mrg_num_, top->P().Eta(), "Y", "top tag eff num y");
                    xBinTopDen = findBin(topTagMisHisto_Mrg_den_, top->P().Pt(),  "X", "top tag eff den x");
                    yBinTopDen = findBin(topTagMisHisto_Mrg_den_, top->P().Eta(), "Y", "top tag eff den y");

                    if (xBinTopNum != -1 and yBinTopNum != -1 and xBinTopDen != -1 and yBinTopDen != -1)
                    {
                        num    = topTagMisHisto_Mrg_num_->GetBinContent(xBinTopNum, yBinTopNum);
                        numUnc = topTagMisHisto_Mrg_num_->GetBinError(xBinTopNum,   yBinTopNum);
                        den    = topTagMisHisto_Mrg_den_->GetBinContent(xBinTopDen, yBinTopDen);
                        denUnc = topTagMisHisto_Mrg_den_->GetBinError(xBinTopDen,   yBinTopDen);
                    }

                    binTopSF = findBin(topMistagSFHisto_Mrg_,     top->P().Pt(), "X", "top tag sf x");
                    if (binTopSF != -1)
                    {
                        sf    = topMistagSFHisto_Mrg_->GetBinContent(binTopSF);
                        sfUnc = topMistagSFHisto_Mrg_->GetBinError(binTopSF);
                    }
                }
            }

            double eff = 1.0, effUnc = 0.0, effUp = 1.0, effDown = 1.0, sfUp = 1.0, sfDown = 1.0;
            if (den != 0.0)
            {
                eff     = num / den;

                // Calculate uncertainty on efficiency ratio in binomial fashion
                // https://root.cern.ch/doc/master/TH1_8cxx_source.html#l03013
                // When num goes to zero, eff goes to zero and effUnc as well
                effUnc  = pow(abs(((1.0 - 2.0*eff)*pow(numUnc, 2.0) + pow(eff, 2.0)*pow(denUnc, 2.0))/pow(den, 2.0)), 0.5);

                effUp   = eff + effUnc;
                effDown = eff - effUnc;

                sfUp    = sf + sfUnc;
                sfDown  = sf - sfUnc;
            }

            if ((isResolved and top->getDiscriminator() > resolvedWP) or
                (isMerged   and top->getDiscriminator() > mergedWP))
            {
                mcTag       *= eff;
                dataTag     *= eff * sf;

                mcTagUp     *= effUp;
                dataTagUp   *= effUp * sfUp;

                mcTagDown   *= effDown;
                dataTagDown *= effDown * sfDown;
            }
            else if (isResolved or isMerged)
            {
                mcNoTag       *= (1.0 - eff);
                dataNoTag     *= (1.0 - eff * sf);

                mcNoTagUp     *= (1.0 - effUp);
                dataNoTagUp   *= (1.0 - effUp * sfUp);

                mcNoTagDown   *= (1.0 - effDown);
                dataNoTagDown *= (1.0 - effDown * sfDown);
            }
        }

        double topTaggerScaleFactor     = (mcNoTag     * mcTag     == 0) ? 1.0 : (dataNoTag     * dataTag    ) / (mcNoTag     * mcTag    );
        double topTaggerScaleFactorUp   = (mcNoTagUp   * mcTagUp   == 0) ? 1.0 : (dataNoTagUp   * dataTagUp  ) / (mcNoTagUp   * mcTagUp  );
        double topTaggerScaleFactorDown = (mcNoTagDown * mcTagDown == 0) ? 1.0 : (dataNoTagDown * dataTagDown) / (mcNoTagDown * mcTagDown);

        tr.registerDerivedVar("topTaggerScaleFactor"     + myVarSuffix_, topTaggerScaleFactor    );
        tr.registerDerivedVar("topTaggerScaleFactorUp"   + myVarSuffix_, topTaggerScaleFactorUp  );
        tr.registerDerivedVar("topTaggerScaleFactorDown" + myVarSuffix_, topTaggerScaleFactorDown);

        //---------------------------------------------------------------------------------------------------------
        // Adding a scale factor for pileup, which comes directly from the nTuples
        // --------------------------------------------------------------------------------------------------------
        const auto& puWeightUnCorr  = tr.getVar<float>("puWeight" );
        const auto& puSysUpUnCorr   = tr.getVar<float>("puSysUp"  );
        const auto& puSysDownUnCorr = tr.getVar<float>("puSysDown");

        const auto& isSignal        = tr.getVar<bool>("isSignal");

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
        double bTagWeightUp = 1.0;
        double bTagWeightDown = 1.0;

        if ( analyzer != "CalculateSFMean" )
        {
            bTagWeight = tr.getVar<double>("bTagSF_EventWeightSimple_Central" +myVarSuffix_);
            bTagWeightUp = tr.getVar<double>("bTagSF_EventWeightSimple_Up" +myVarSuffix_);
            bTagWeightDown = tr.getVar<double>("bTagSF_EventWeightSimple_Down" +myVarSuffix_);
        }

        double CommonWeight      = Weight * FinalLumi * topPtScaleFactor;
        double CommonWeightQCDCR = totNonIsoMuonSF * prefiringScaleFactor * puWeightCorr;
        double CommonWeight0l    = jetTrigSF * bTagWeight * prefiringScaleFactor * puWeightCorr * topTaggerScaleFactor;
        double CommonWeight1l    = totGoodElectronSF * totGoodMuonSF * bTagWeight * prefiringScaleFactor * puWeightCorr;
        double CommonWeight2l    = totGoodElectronSF * totGoodMuonSF * bTagWeight * prefiringScaleFactor * puWeightCorr;

        double totalEventWeight_0l         = CommonWeight * CommonWeight0l;
        double totalEventWeight_0l_TtgUp   = CommonWeight * topTaggerScaleFactorUp   * jetTrigSF      * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_TtgDown = CommonWeight * topTaggerScaleFactorDown * jetTrigSF      * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_JetUp   = CommonWeight * topTaggerScaleFactor     * jetTrigSF_Up   * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_JetDown = CommonWeight * topTaggerScaleFactor     * jetTrigSF_Down * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_BtgUp   = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeightUp   * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_BtgDown = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeightDown * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_0l_PrfUp   = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeight     * prefiringScaleFactorUp   * puWeightCorr;
        double totalEventWeight_0l_PrfDown = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeight     * prefiringScaleFactorDown * puWeightCorr;
        double totalEventWeight_0l_PUup    = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeight     * prefiringScaleFactor     * puSysUpCorr;
        double totalEventWeight_0l_PUdown  = CommonWeight * topTaggerScaleFactor     * jetTrigSF      * bTagWeight     * prefiringScaleFactor     * puSysDownCorr;
        double totalEventWeight_0l_SclUp   = CommonWeight * CommonWeight0l * scaleWeightUpperBound_corr;
        double totalEventWeight_0l_SclDown = CommonWeight * CommonWeight0l * scaleWeightLowerBound_corr;
        double totalEventWeight_0l_PDFup   = CommonWeight * CommonWeight0l * NNPDF_from_median_up_corr;
        double totalEventWeight_0l_PDFdown = CommonWeight * CommonWeight0l * NNPDF_from_median_down_corr;
        double totalEventWeight_0l_ISRup   = CommonWeight * CommonWeight0l * PSweight_ISRUp;
        double totalEventWeight_0l_ISRdown = CommonWeight * CommonWeight0l * PSweight_ISRDown;
        double totalEventWeight_0l_FSRup   = CommonWeight * CommonWeight0l * PSweight_FSRUp;
        double totalEventWeight_0l_FSRdown = CommonWeight * CommonWeight0l * PSweight_FSRDown;

        double totalEventWeight_1l         = CommonWeight * CommonWeight1l;
        double totalEventWeight_1l_BtgUp   = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeightUp   * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_1l_BtgDown = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeightDown * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_1l_LepUp   = CommonWeight * totGoodElectronSF_Up   * totGoodMuonSF_Up   * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_1l_LepDown = CommonWeight * totGoodElectronSF_Down * totGoodMuonSF_Down * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_1l_PUup    = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysUpCorr;
        double totalEventWeight_1l_PUdown  = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysDownCorr;
        double totalEventWeight_1l_PrfUp   = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactorUp   * puWeightCorr;
        double totalEventWeight_1l_PrfDown = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactorDown * puWeightCorr;
        double totalEventWeight_1l_SclUp   = CommonWeight * CommonWeight1l * scaleWeightUpperBound_corr;
        double totalEventWeight_1l_SclDown = CommonWeight * CommonWeight1l * scaleWeightLowerBound_corr;
        double totalEventWeight_1l_PDFup   = CommonWeight * CommonWeight1l * NNPDF_from_median_up_corr;
        double totalEventWeight_1l_PDFdown = CommonWeight * CommonWeight1l * NNPDF_from_median_down_corr;
        double totalEventWeight_1l_ISRup   = CommonWeight * CommonWeight1l * PSweight_ISRUp;
        double totalEventWeight_1l_ISRdown = CommonWeight * CommonWeight1l * PSweight_ISRDown;
        double totalEventWeight_1l_FSRup   = CommonWeight * CommonWeight1l * PSweight_FSRUp;
        double totalEventWeight_1l_FSRdown = CommonWeight * CommonWeight1l * PSweight_FSRDown;

        double totalEventWeight_2l         = CommonWeight * CommonWeight2l;
        double totalEventWeight_2l_BtgUp   = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeightUp   * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_2l_BtgDown = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeightDown * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_2l_LepUp   = CommonWeight * totGoodElectronSF_Up   * totGoodMuonSF_Up   * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_2l_LepDown = CommonWeight * totGoodElectronSF_Down * totGoodMuonSF_Down * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_2l_PUup    = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysUpCorr;
        double totalEventWeight_2l_PUdown  = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysDownCorr;
        double totalEventWeight_2l_PrfUp   = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactorUp   * puWeightCorr;
        double totalEventWeight_2l_PrfDown = CommonWeight * totGoodElectronSF      * totGoodMuonSF      * bTagWeight     * prefiringScaleFactorDown * puWeightCorr;
        double totalEventWeight_2l_SclUp   = CommonWeight * CommonWeight2l * scaleWeightUpperBound_corr;
        double totalEventWeight_2l_SclDown = CommonWeight * CommonWeight2l * scaleWeightLowerBound_corr;
        double totalEventWeight_2l_PDFup   = CommonWeight * CommonWeight2l * NNPDF_from_median_up_corr;
        double totalEventWeight_2l_PDFdown = CommonWeight * CommonWeight2l * NNPDF_from_median_down_corr;
        double totalEventWeight_2l_ISRup   = CommonWeight * CommonWeight2l * PSweight_ISRUp;
        double totalEventWeight_2l_ISRdown = CommonWeight * CommonWeight2l * PSweight_ISRDown;
        double totalEventWeight_2l_FSRup   = CommonWeight * CommonWeight2l * PSweight_FSRUp;
        double totalEventWeight_2l_FSRdown = CommonWeight * CommonWeight2l * PSweight_FSRDown;

        double totalEventWeight_QCDCR          = CommonWeight * CommonWeightQCDCR;
        double totalEventWeight_QCDCR_NimUp    = CommonWeight * totNonIsoMuonSF_Up   * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_QCDCR_NimDown  = CommonWeight * totNonIsoMuonSF_Down * bTagWeight     * prefiringScaleFactor     * puWeightCorr;
        double totalEventWeight_QCDCR_PrfUp    = CommonWeight * totNonIsoMuonSF      * bTagWeight     * prefiringScaleFactorUp   * puWeightCorr;
        double totalEventWeight_QCDCR_PrfDown  = CommonWeight * totNonIsoMuonSF      * bTagWeight     * prefiringScaleFactorDown * puWeightCorr;
        double totalEventWeight_QCDCR_PUup     = CommonWeight * totNonIsoMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysUpCorr;
        double totalEventWeight_QCDCR_PUdown   = CommonWeight * totNonIsoMuonSF      * bTagWeight     * prefiringScaleFactor     * puSysDownCorr;
        double totalEventWeight_QCDCR_SclUp    = CommonWeight * CommonWeightQCDCR * scaleWeightUpperBound_corr;
        double totalEventWeight_QCDCR_SclDown  = CommonWeight * CommonWeightQCDCR * scaleWeightLowerBound_corr;
        double totalEventWeight_QCDCR_PDFup    = CommonWeight * CommonWeightQCDCR * NNPDF_from_median_up_corr;
        double totalEventWeight_QCDCR_PDFdown  = CommonWeight * CommonWeightQCDCR * NNPDF_from_median_down_corr;
        double totalEventWeight_QCDCR_ISRup    = CommonWeight * CommonWeightQCDCR * PSweight_ISRUp;
        double totalEventWeight_QCDCR_ISRdown  = CommonWeight * CommonWeightQCDCR * PSweight_ISRDown;
        double totalEventWeight_QCDCR_FSRup    = CommonWeight * CommonWeightQCDCR * PSweight_FSRUp;
        double totalEventWeight_QCDCR_FSRdown  = CommonWeight * CommonWeightQCDCR * PSweight_FSRDown;

        tr.registerDerivedVar("TotalWeight_QCDCR"          + myVarSuffix_, totalEventWeight_QCDCR);
        tr.registerDerivedVar("TotalWeight_QCDCR_NimUp"    + myVarSuffix_, totalEventWeight_QCDCR_NimUp);
        tr.registerDerivedVar("TotalWeight_QCDCR_NimDown"  + myVarSuffix_, totalEventWeight_QCDCR_NimDown);
        tr.registerDerivedVar("TotalWeight_QCDCR_PUup"     + myVarSuffix_, totalEventWeight_QCDCR_PUup);
        tr.registerDerivedVar("TotalWeight_QCDCR_PUdown"   + myVarSuffix_, totalEventWeight_QCDCR_PUdown);
        tr.registerDerivedVar("TotalWeight_QCDCR_PrfUp"    + myVarSuffix_, totalEventWeight_QCDCR_PrfUp);
        tr.registerDerivedVar("TotalWeight_QCDCR_PrfDown"  + myVarSuffix_, totalEventWeight_QCDCR_PrfDown);
        tr.registerDerivedVar("TotalWeight_QCDCR_SclUp"    + myVarSuffix_, totalEventWeight_QCDCR_SclUp);
        tr.registerDerivedVar("TotalWeight_QCDCR_SclDown"  + myVarSuffix_, totalEventWeight_QCDCR_SclDown);
        tr.registerDerivedVar("TotalWeight_QCDCR_PDFup"    + myVarSuffix_, totalEventWeight_QCDCR_PDFup);
        tr.registerDerivedVar("TotalWeight_QCDCR_PDFdown"  + myVarSuffix_, totalEventWeight_QCDCR_PDFdown);
        tr.registerDerivedVar("TotalWeight_QCDCR_ISRup"    + myVarSuffix_, totalEventWeight_QCDCR_ISRup);
        tr.registerDerivedVar("TotalWeight_QCDCR_ISRdown"  + myVarSuffix_, totalEventWeight_QCDCR_ISRdown);
        tr.registerDerivedVar("TotalWeight_QCDCR_FSRup"    + myVarSuffix_, totalEventWeight_QCDCR_FSRup);
        tr.registerDerivedVar("TotalWeight_QCDCR_FSRdown"  + myVarSuffix_, totalEventWeight_QCDCR_FSRdown);

        tr.registerDerivedVar("TotalWeight_0l"           + myVarSuffix_, totalEventWeight_0l);
        tr.registerDerivedVar("TotalWeight_0l_BtgUp"     + myVarSuffix_, totalEventWeight_0l_BtgUp);
        tr.registerDerivedVar("TotalWeight_0l_BtgDown"   + myVarSuffix_, totalEventWeight_0l_BtgDown);
        tr.registerDerivedVar("TotalWeight_0l_TtgUp"     + myVarSuffix_, totalEventWeight_0l_TtgUp);
        tr.registerDerivedVar("TotalWeight_0l_TtgDown"   + myVarSuffix_, totalEventWeight_0l_TtgDown);
        tr.registerDerivedVar("TotalWeight_0l_JetUp"     + myVarSuffix_, totalEventWeight_0l_JetUp);
        tr.registerDerivedVar("TotalWeight_0l_JetDown"   + myVarSuffix_, totalEventWeight_0l_JetDown);
        tr.registerDerivedVar("TotalWeight_0l_PUup"      + myVarSuffix_, totalEventWeight_0l_PUup);
        tr.registerDerivedVar("TotalWeight_0l_PUdown"    + myVarSuffix_, totalEventWeight_0l_PUdown);
        tr.registerDerivedVar("TotalWeight_0l_PrfUp"     + myVarSuffix_, totalEventWeight_0l_PrfUp);
        tr.registerDerivedVar("TotalWeight_0l_PrfDown"   + myVarSuffix_, totalEventWeight_0l_PrfDown);
        tr.registerDerivedVar("TotalWeight_0l_SclUp"     + myVarSuffix_, totalEventWeight_0l_SclUp);
        tr.registerDerivedVar("TotalWeight_0l_SclDown"   + myVarSuffix_, totalEventWeight_0l_SclDown);
        tr.registerDerivedVar("TotalWeight_0l_PDFup"     + myVarSuffix_, totalEventWeight_0l_PDFup);
        tr.registerDerivedVar("TotalWeight_0l_PDFdown"   + myVarSuffix_, totalEventWeight_0l_PDFdown);
        tr.registerDerivedVar("TotalWeight_0l_ISRup"     + myVarSuffix_, totalEventWeight_0l_ISRup);
        tr.registerDerivedVar("TotalWeight_0l_ISRdown"   + myVarSuffix_, totalEventWeight_0l_ISRdown);
        tr.registerDerivedVar("TotalWeight_0l_FSRup"     + myVarSuffix_, totalEventWeight_0l_FSRup);
        tr.registerDerivedVar("TotalWeight_0l_FSRdown"   + myVarSuffix_, totalEventWeight_0l_FSRdown);

        tr.registerDerivedVar("TotalWeight_1l"           + myVarSuffix_, totalEventWeight_1l);
        tr.registerDerivedVar("TotalWeight_1l_BtgUp"     + myVarSuffix_, totalEventWeight_1l_BtgUp);
        tr.registerDerivedVar("TotalWeight_1l_BtgDown"   + myVarSuffix_, totalEventWeight_1l_BtgDown);
        tr.registerDerivedVar("TotalWeight_1l_LepUp"     + myVarSuffix_, totalEventWeight_1l_LepUp);
        tr.registerDerivedVar("TotalWeight_1l_LepDown"   + myVarSuffix_, totalEventWeight_1l_LepDown);
        tr.registerDerivedVar("TotalWeight_1l_PUup"      + myVarSuffix_, totalEventWeight_1l_PUup);
        tr.registerDerivedVar("TotalWeight_1l_PUdown"    + myVarSuffix_, totalEventWeight_1l_PUdown);
        tr.registerDerivedVar("TotalWeight_1l_PrfUp"     + myVarSuffix_, totalEventWeight_1l_PrfUp);
        tr.registerDerivedVar("TotalWeight_1l_PrfDown"   + myVarSuffix_, totalEventWeight_1l_PrfDown);
        tr.registerDerivedVar("TotalWeight_1l_SclUp"     + myVarSuffix_, totalEventWeight_1l_SclUp);
        tr.registerDerivedVar("TotalWeight_1l_SclDown"   + myVarSuffix_, totalEventWeight_1l_SclDown);
        tr.registerDerivedVar("TotalWeight_1l_PDFup"     + myVarSuffix_, totalEventWeight_1l_PDFup);
        tr.registerDerivedVar("TotalWeight_1l_PDFdown"   + myVarSuffix_, totalEventWeight_1l_PDFdown);
        tr.registerDerivedVar("TotalWeight_1l_ISRup"     + myVarSuffix_, totalEventWeight_1l_ISRup);
        tr.registerDerivedVar("TotalWeight_1l_ISRdown"   + myVarSuffix_, totalEventWeight_1l_ISRdown);
        tr.registerDerivedVar("TotalWeight_1l_FSRup"     + myVarSuffix_, totalEventWeight_1l_FSRup);
        tr.registerDerivedVar("TotalWeight_1l_FSRdown"   + myVarSuffix_, totalEventWeight_1l_FSRdown);

        tr.registerDerivedVar("TotalWeight_2l"           + myVarSuffix_, totalEventWeight_2l);
        tr.registerDerivedVar("TotalWeight_2l_BtgUp"     + myVarSuffix_, totalEventWeight_2l_BtgUp);
        tr.registerDerivedVar("TotalWeight_2l_BtgDown"   + myVarSuffix_, totalEventWeight_2l_BtgDown);
        tr.registerDerivedVar("TotalWeight_2l_LepUp"     + myVarSuffix_, totalEventWeight_2l_LepUp);
        tr.registerDerivedVar("TotalWeight_2l_LepDown"   + myVarSuffix_, totalEventWeight_2l_LepDown);
        tr.registerDerivedVar("TotalWeight_2l_PUup"      + myVarSuffix_, totalEventWeight_2l_PUup);
        tr.registerDerivedVar("TotalWeight_2l_PUdown"    + myVarSuffix_, totalEventWeight_2l_PUdown);
        tr.registerDerivedVar("TotalWeight_2l_PrfUp"     + myVarSuffix_, totalEventWeight_2l_PrfUp);
        tr.registerDerivedVar("TotalWeight_2l_PrfDown"   + myVarSuffix_, totalEventWeight_2l_PrfDown);
        tr.registerDerivedVar("TotalWeight_2l_SclUp"     + myVarSuffix_, totalEventWeight_2l_SclUp);
        tr.registerDerivedVar("TotalWeight_2l_SclDown"   + myVarSuffix_, totalEventWeight_2l_SclDown);
        tr.registerDerivedVar("TotalWeight_2l_PDFup"     + myVarSuffix_, totalEventWeight_2l_PDFup);
        tr.registerDerivedVar("TotalWeight_2l_PDFdown"   + myVarSuffix_, totalEventWeight_2l_PDFdown);
        tr.registerDerivedVar("TotalWeight_2l_ISRup"     + myVarSuffix_, totalEventWeight_2l_ISRup);
        tr.registerDerivedVar("TotalWeight_2l_ISRdown"   + myVarSuffix_, totalEventWeight_2l_ISRdown);
        tr.registerDerivedVar("TotalWeight_2l_FSRup"     + myVarSuffix_, totalEventWeight_2l_FSRup);
        tr.registerDerivedVar("TotalWeight_2l_FSRdown"   + myVarSuffix_, totalEventWeight_2l_FSRdown);

        if(printGetBinError_) firstPrint_ = false;
    }

public:
    ScaleFactors( const std::string& runYear, const std::string& leptonFileName, const std::string& hadronicFileName, const std::string& toptaggerFileName, const std::string& meanFileName, const std::string& filetag, const std::string& myVarSuffix = "" )
        : myVarSuffix_(myVarSuffix), firstPrint_(true), printMeanError_(false), printGetBinError_(false)
    {
        std::cout<<"Setting up ScaleFactors"<<std::endl;

        // If we run on a skim data set, there should be a "_skim" suffix
        // at the end of the name, so scrape it off before trying
        // to get anything that requires the filetag
        std::string scrapeOff    = "_skim";
        std::string filetagClean = filetag;
        const size_t position = filetag.find(scrapeOff);
        if(position != std::string::npos)
        {
            filetagClean.erase(position, scrapeOff.length());
        }

        // Force histograms to reside in memory rather than a TFile buffer on disk, allowing the TFile to be closed safely
        TH1::AddDirectory(false);

        // Getting Leptonic and Hadronic scale factor histograms
        TFile leptonic_SFRootFile( leptonFileName.c_str()   );
        TFile hadronic_SFRootFile( hadronicFileName.c_str() );
        TFile toptagger_SFRootFile( toptaggerFileName.c_str() );

        TString topTagSFHistoName_Res              = runYear + "_TagRateSF_vs_topPt_Resolved";
        TString topTagSFHistoName_Mrg              = runYear + "_TagRateSF_vs_topPt_Merged";
        TString topMistagSFHistoName_Res           = runYear + "_MisTagSF_vs_topPt_Resolved";
        TString topMistagSFHistoName_Mrg           = runYear + "_MisTagSF_vs_topPt_Merged";
        TString topTagEffHistoName_Mrg_den         = "d_eff_mrg_" + filetagClean;
        TString topTagEffHistoName_Res_den         = "d_eff_res_" + filetagClean;
        TString topTagEffHistoName_Mrg_num         = "n_eff_mrg_" + filetagClean;
        TString topTagEffHistoName_Res_num         = "n_eff_res_" + filetagClean;
        TString topTagMisHistoName_Mrg_den         = "d_mis_mrg_" + filetagClean;
        TString topTagMisHistoName_Res_den         = "d_mis_res_" + filetagClean;
        TString topTagMisHistoName_Mrg_num         = "n_mis_mrg_" + filetagClean;
        TString topTagMisHistoName_Res_num         = "n_mis_res_" + filetagClean;
        TString eleSFHistoTightName                = "EGamma_SF2D_" + runYear + "_UL_ID";
        TString eleSFHistoRecoName                 = "EGamma_SF2D_" + runYear + "_UL_RECO";
        TString eleSFHistoTrigName                 = runYear + "_el_pt40_trig_ge5jetCut_wLepPtLepEtaBin_TriggerSF";
        TString muSFHistoMediumName                = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_" + runYear + "_UL_ID";
        TString muSFHistoIsoName                   = "NUM_TightRelIso_DEN_MediumID_abseta_pt_" + runYear + "_UL_ISO";
        TString muSFHistoTrigName                  = runYear + "_mu_pt40_trig_ge5jetCut_wLepPtLepEtaBin_TriggerSF";
        TString nimuSFHistoTrigName                = ""; //just for calculating non iso muon scale factors
        TString jetSFHistoTrigName_1bCut           = runYear + "_jet_trig_1bjetCut_wJetHt6thJetPtBin_TriggerSF";
        TString jetSFHistoTrigName_ge2bCut         = runYear + "_jet_trig_ge2bjetCut_wJetHt6thJetPtBin_TriggerSF"; 

        getHisto(leptonic_SFRootFile, eleSFHistoTight_,            eleSFHistoTightName       );
        getHisto(leptonic_SFRootFile, eleSFHistoReco_,             eleSFHistoRecoName        );
        getHisto(leptonic_SFRootFile, eleSFHistoTrig_,             eleSFHistoTrigName        );
        getHisto(leptonic_SFRootFile, muSFHistoMedium_,            muSFHistoMediumName       );
        getHisto(leptonic_SFRootFile, muSFHistoIso_,               muSFHistoIsoName          );
        getHisto(leptonic_SFRootFile, muSFHistoTrig_,              muSFHistoTrigName         );
        getHisto(leptonic_SFRootFile, nimuSFHistoTrig_,            nimuSFHistoTrigName       );
        getHisto(hadronic_SFRootFile, jetSFHistoTrigName_1bCut_,   jetSFHistoTrigName_1bCut  );
        getHisto(hadronic_SFRootFile, jetSFHistoTrigName_ge2bCut_, jetSFHistoTrigName_ge2bCut);        
        getHisto(toptagger_SFRootFile, topTagSFHisto_Res_,              topTagSFHistoName_Res             );
        getHisto(toptagger_SFRootFile, topTagSFHisto_Mrg_,              topTagSFHistoName_Mrg             );
        getHisto(toptagger_SFRootFile, topMistagSFHisto_Res_,           topMistagSFHistoName_Res          );
        getHisto(toptagger_SFRootFile, topMistagSFHisto_Mrg_,           topMistagSFHistoName_Mrg          );
        getHisto(toptagger_SFRootFile, topTagEffHisto_Mrg_den_,         topTagEffHistoName_Mrg_den        );
        getHisto(toptagger_SFRootFile, topTagEffHisto_Res_den_,         topTagEffHistoName_Res_den        );
        getHisto(toptagger_SFRootFile, topTagMisHisto_Mrg_den_,         topTagMisHistoName_Mrg_den        );
        getHisto(toptagger_SFRootFile, topTagMisHisto_Res_den_,         topTagMisHistoName_Res_den        );
        getHisto(toptagger_SFRootFile, topTagEffHisto_Mrg_num_,         topTagEffHistoName_Mrg_num        );
        getHisto(toptagger_SFRootFile, topTagEffHisto_Res_num_,         topTagEffHistoName_Res_num        );
        getHisto(toptagger_SFRootFile, topTagMisHisto_Mrg_num_,         topTagMisHistoName_Mrg_num        );
        getHisto(toptagger_SFRootFile, topTagMisHisto_Res_num_,         topTagMisHistoName_Res_num        );

        leptonic_SFRootFile.Close();
        hadronic_SFRootFile.Close();
        toptagger_SFRootFile.Close();

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

        if (!lostCauseEvent or !fastMode or tr.isFirstEvent())
            scaleFactors(tr);
    }
};

#endif
