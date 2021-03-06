#ifndef RUNFISHER_H
#define RUNFISHER_H

// includes for the event shapes
#include "Framework/Framework/include/EventShapeVariables.h"
#include "Framework/Framework/src/get_cmframe_jets.c"
#include "Framework/Framework/include/bdt_350to650_fwm10_jmtev_top6.h"
#include "Framework/Framework/src/fisher_350to650_fwm10_jmtev_top6.c"
#include "Framework/Framework/src/fisher_350to650_fwm6_jmtev_top6_gt_v2.c"
#include "Framework/Framework/src/fisher_350to650_fwm6_jmtev_top6_gt_v3pt30.c"
#include "Framework/Framework/src/fisher_0lepton_v1.c"

//#include "Framework/BackgroundMVA/test/fisherLoader/weights/TMVAClassification_FisherG.class.C"

class RunFisher
{
private:
    std::string myVarSuffix_;
    std::string fisherVersion_;
    std::vector<std::string> inputVarNames_top6_fwm6_;
    std::vector<std::string> inputVarNames_top6_fwm10_;
    std::vector<double> inputVals_top6_fwm6_;
    std::vector<double> inputVals_top6_fwm10_;
    std::shared_ptr<ReadBDT_350to650_fwm10_jmtev_top6>              eventshapeBDT_;
    std::shared_ptr<ReadFisher_350to650_fwm10_jmtev_top6>           read_fisher_350to650_fwm10_jmtev_top6_;
    std::shared_ptr<ReadFisherG_350to650_fwm6_jmtev_top6_gt_v2>     read_fisher_350to650_fwm6_jmtev_top6_gt_v2_;
    std::shared_ptr<ReadFisherG_350to650_fwm6_jmtev_top6_gt_v3pt30> read_fisher_350to650_fwm6_jmtev_top6_gt_v3pt30_;
    std::shared_ptr<ReadFisherG_0lepton_v1>                         read_fisher_0lepton_v1_;

    //std::shared_ptr<ReadFisherG>                                    read_fisher_test_;

    void prepareFWMVecs(const std::vector<std::string>& names, std::vector<double>& values)
    {
        for ( unsigned int i=0; i < names.size() ; i++ ) 
        {
            values.push_back( 0.5 ) ; //--- load vector with dummy values.
        }
    }

    void resetFWMVecs(std::vector<double>& v, EventShapeVariables& esv)
    {
        int index = -1;
        TVectorD eigen_vals_norm_top6 = esv.getEigenValues();
        for(int i = 0; i < v.size(); i++)
        {
            if(i < v.size() - 3)
            {
                v[i] = esv.getFWmoment(i+2);
            }
            else
            {
                index++;
                v[i] = eigen_vals_norm_top6[index];
            }
        }
    }

    void setUpFWM()
    {
        // FWM 2-6
        inputVarNames_top6_fwm6_ = {"fwm2_top6","fwm3_top6","fwm4_top6","fwm5_top6","fwm6_top6",
                                    "jmt_ev0_top6","jmt_ev1_top6","jmt_ev2_top6"};
        
        prepareFWMVecs(inputVarNames_top6_fwm6_, inputVals_top6_fwm6_);

        // FWM 2-10
        inputVarNames_top6_fwm10_  = {"fwm2_top6","fwm3_top6","fwm4_top6","fwm5_top6","fwm6_top6","fwm7_top6","fwm8_top6","fwm9_top6","fwm10_top6",
                                      "jmt_ev0_top6","jmt_ev1_top6","jmt_ev2_top6"};

        prepareFWMVecs(inputVarNames_top6_fwm10_, inputVals_top6_fwm10_);
    }

    void runFisher(NTupleReader& tr)
    {
        const auto& Jets       = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& NJets_pt30 = tr.getVar<int>("NJets_pt30"+myVarSuffix_);
        const auto& NJets_pt45 = tr.getVar<int>("NJets_pt45"+myVarSuffix_);

        std::vector<math::RThetaPhiVector> cm_frame_jets;
        get_cmframe_jets( &Jets, cm_frame_jets, 6 );
        EventShapeVariables esv_top6( cm_frame_jets );
        resetFWMVecs(inputVals_top6_fwm10_, esv_top6);
        resetFWMVecs(inputVals_top6_fwm6_, esv_top6);
        double eventshape_bdt_val = eventshapeBDT_->GetMvaValue( inputVals_top6_fwm10_ );        

        double fisher_val;
        if( fisherVersion_ == "v1" )
        {            
            fisher_val = read_fisher_350to650_fwm10_jmtev_top6_->GetMvaValue( inputVals_top6_fwm10_ );
        }
        else if ( fisherVersion_ == "v2" )
        {
            fisher_val = read_fisher_350to650_fwm6_jmtev_top6_gt_v2_->GetMvaValue( inputVals_top6_fwm6_ );
        }
        else if ( fisherVersion_ == "v3" )
        {
            fisher_val = read_fisher_350to650_fwm6_jmtev_top6_gt_v3pt30_->GetMvaValue( inputVals_top6_fwm6_ );
            // Apply the fisher correction in bins of njets
            std::map<int, double> fisher_shift = {
                {6, -0.0001},
                {7, -0.0003},
                {8,  0.0008},
                {9,  0.0018},
                {10, 0.0030},
                {11, 0.0038},
                {12, 0.0043},
                {13, 0.0086},
                {14,-0.0026},
                {15, 0.0207},
            };

            std::map<int, double> fisher_width_shift = {
                {6,  0.99248},
                {7,  1.00053},
                {8,  1.00014},
                {9,  0.99957},
                {10, 0.99844},
                {11, 1.00922},
                {12, 1.00409},
                {13, 1.02699},
                {14, 1.01599},
                {15, 1.00572},
            };
            if (NJets_pt30 >= 6 && NJets_pt30 <= 16)
            {
                //fisher_val = fisher_val + fisher_shift[NJets_pt30];
                fisher_val = fisher_width_shift[NJets_pt30]*(fisher_val + fisher_shift[NJets_pt30]);
            }
        }
        else if ( fisherVersion_ == "0lepton_v1")
        {
            fisher_val = read_fisher_0lepton_v1_->GetMvaValue( inputVals_top6_fwm6_ );
            // Apply the fisher correction in bins of njets
            std::map<int, double> fisher_shift = {
                {6,  -0.0001644},
                {7,  -0.0002359},
                {8,   0.0004749},
                {9,   0.0015723},
                {10,  0.0002881},
                {11,  0.0027997},
                {12, -0.0003774},
                {13, -0.0136242},
                {14, -0.0368304},
                {15,  0.0209468},
            };
            if (NJets_pt30 >= 6 && NJets_pt30 <= 16)
            {
                fisher_val = fisher_val + fisher_shift[NJets_pt30];
            }
        }
        //else if ( fisherVersion_ == "test")
        //{
        //    fisher_val = read_fisher_test_->GetMvaValue( inputVals_top6_fwm6_ );
        //}

        bool bdt_bin1 = eventshape_bdt_val > -1.00 && eventshape_bdt_val <= -0.04;
        bool bdt_bin2 = eventshape_bdt_val > -0.04 && eventshape_bdt_val <=  0.00;
        bool bdt_bin3 = eventshape_bdt_val >  0.00 && eventshape_bdt_val <=  0.04;
        bool bdt_bin4 = eventshape_bdt_val >  0.04 && eventshape_bdt_val <=  1.00;

        bool fisher_bin1 = fisher_val > -1.000 && fisher_val <= -0.035;
        bool fisher_bin2 = fisher_val > -0.035 && fisher_val <=  0.030;
        bool fisher_bin3 = fisher_val >  0.030 && fisher_val <=  0.095;
        bool fisher_bin4 = fisher_val >  0.095 && fisher_val <=  1.000;

        if (fisherVersion_ == "v3")
        {
            fisher_bin1 = fisher_val > -1.000 && fisher_val <= -0.015;
            fisher_bin2 = fisher_val > -0.015 && fisher_val <=  0.020;
            fisher_bin3 = fisher_val >  0.020 && fisher_val <=  0.060;
            fisher_bin4 = fisher_val >  0.060 && fisher_val <=  1.000;            
        }
        else if (fisherVersion_ == "0lepton_v1")
        {
            fisher_bin1 = fisher_val > -1.0000 && fisher_val <= -0.1225;
            fisher_bin2 = fisher_val > -0.1225 && fisher_val <= -0.0725;
            fisher_bin3 = fisher_val > -0.0725 && fisher_val <= -0.0175;
            fisher_bin4 = fisher_val > -0.0175 && fisher_val <=  1.0000;            
        }

        // Register Variables
        tr.registerDerivedVar("eventshape_bdt_val"+myVarSuffix_, eventshape_bdt_val);
        tr.registerDerivedVar("bdt_bin1"+myVarSuffix_, bdt_bin1);
        tr.registerDerivedVar("bdt_bin2"+myVarSuffix_, bdt_bin2);
        tr.registerDerivedVar("bdt_bin3"+myVarSuffix_, bdt_bin3);
        tr.registerDerivedVar("bdt_bin4"+myVarSuffix_, bdt_bin4);
        tr.registerDerivedVar("fisher_val"+myVarSuffix_, fisher_val);
        tr.registerDerivedVar("fisher_bin1"+myVarSuffix_, fisher_bin1);
        tr.registerDerivedVar("fisher_bin2"+myVarSuffix_, fisher_bin2);
        tr.registerDerivedVar("fisher_bin3"+myVarSuffix_, fisher_bin3);
        tr.registerDerivedVar("fisher_bin4"+myVarSuffix_, fisher_bin4);
    }

public:
    RunFisher(std::string fisherVersion = "v3", std::string myVarSuffix = "") 
        : fisherVersion_(fisherVersion)
        , myVarSuffix_(myVarSuffix)
        , eventshapeBDT_(nullptr)
        , read_fisher_350to650_fwm10_jmtev_top6_(nullptr)
        , read_fisher_350to650_fwm6_jmtev_top6_gt_v2_(nullptr)
        , read_fisher_350to650_fwm6_jmtev_top6_gt_v3pt30_(nullptr)
        , read_fisher_0lepton_v1_(nullptr)
        //, read_fisher_test_(nullptr)
    {
        std::cout<<"Setting up RunFisher"<<std::endl;
        setUpFWM();
        eventshapeBDT_                                  = std::make_shared<ReadBDT_350to650_fwm10_jmtev_top6>( inputVarNames_top6_fwm10_ );
        read_fisher_350to650_fwm10_jmtev_top6_          = std::make_shared<ReadFisher_350to650_fwm10_jmtev_top6>( inputVarNames_top6_fwm10_ );
        read_fisher_350to650_fwm6_jmtev_top6_gt_v2_     = std::make_shared<ReadFisherG_350to650_fwm6_jmtev_top6_gt_v2>( inputVarNames_top6_fwm6_ );
        read_fisher_350to650_fwm6_jmtev_top6_gt_v3pt30_ = std::make_shared<ReadFisherG_350to650_fwm6_jmtev_top6_gt_v3pt30>( inputVarNames_top6_fwm6_ );
        read_fisher_0lepton_v1_                         = std::make_shared<ReadFisherG_0lepton_v1>( inputVarNames_top6_fwm6_ );
        //read_fisher_test_                               = std::make_shared<ReadFisherG>( inputVarNames_top6_fwm6_ );
        std::cout<<"Using Fisher version: "+fisherVersion<<std::endl;
    }
    
    void operator()(NTupleReader& tr)
    {
        runFisher(tr);
    }
};

#endif
