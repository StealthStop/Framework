#ifndef DEEPEVENTSHAPE_H
#define DEEPEVENTSHAPE_H

#include "tensorflow/c/c_api.h"
#include "TopTagger/CfgParser/include/Context.hh"
#include "TopTagger/CfgParser/include/CfgDocument.hh"

#include "cstdlib"
#include "cstdio"
#include "cstring"

class EventShapeCalculator
{
private:
    float* basePtr_;
    int len_;

    int fwm2_top6_, fwm3_top6_, fwm4_top6_, fwm5_top6_, fwm6_top6_, fwm7_top6_, fwm8_top6_, fwm9_top6_, fwm10_top6_, jmt_ev0_top6_, jmt_ev1_top6_, jmt_ev2_top6_;
    int NGoodJets_double_;
    int Jet_pt_1_, Jet_pt_2_, Jet_pt_3_, Jet_pt_4_, Jet_pt_5_, Jet_pt_6_, Jet_pt_7_, Jet_pt_8_, Jet_pt_9_, Jet_pt_10_;
    int Jet_pt_11_, Jet_pt_12_, Jet_pt_13_, Jet_pt_14_, Jet_pt_15_, Jet_pt_16_, Jet_pt_17_, Jet_pt_18_, Jet_pt_19_, Jet_pt_20_;
    int Jet_eta_1_, Jet_eta_2_, Jet_eta_3_, Jet_eta_4_, Jet_eta_5_, Jet_eta_6_, Jet_eta_7_, Jet_eta_8_, Jet_eta_9_, Jet_eta_10_;
    int Jet_eta_11_, Jet_eta_12_, Jet_eta_13_, Jet_eta_14_, Jet_eta_15_, Jet_eta_16_, Jet_eta_17_, Jet_eta_18_, Jet_eta_19_, Jet_eta_20_;
    int Jet_phi_1_, Jet_phi_2_, Jet_phi_3_, Jet_phi_4_, Jet_phi_5_, Jet_phi_6_, Jet_phi_7_, Jet_phi_8_, Jet_phi_9_, Jet_phi_10_;
    int Jet_phi_11_, Jet_phi_12_, Jet_phi_13_, Jet_phi_14_, Jet_phi_15_, Jet_phi_16_, Jet_phi_17_, Jet_phi_18_, Jet_phi_19_, Jet_phi_20_;
    int Jet_m_1_, Jet_m_2_, Jet_m_3_, Jet_m_4_, Jet_m_5_, Jet_m_6_, Jet_m_7_, Jet_m_8_, Jet_m_9_, Jet_m_10_;
    int Jet_m_11_, Jet_m_12_, Jet_m_13_, Jet_m_14_, Jet_m_15_, Jet_m_16_, Jet_m_17_, Jet_m_18_, Jet_m_19_, Jet_m_20_;
    int GoodLeptons_pt_1_, GoodLeptons_eta_1_, GoodLeptons_phi_1_, GoodLeptons_m_1_;
    int BestComboAvgMass_;

public:
    EventShapeCalculator()
    {
        fwm2_top6_ = fwm3_top6_ = fwm4_top6_ = fwm5_top6_ = fwm6_top6_ = fwm7_top6_ = fwm8_top6_ = fwm9_top6_ = fwm10_top6_ = jmt_ev0_top6_ = jmt_ev1_top6_ = jmt_ev2_top6_ = -1; 
        NGoodJets_double_ = -1;
        Jet_pt_1_ = Jet_pt_2_ = Jet_pt_3_ = Jet_pt_4_ = Jet_pt_5_ = Jet_pt_6_ = Jet_pt_7_ = Jet_pt_8_ = Jet_pt_9_ = Jet_pt_10_ = -1;
        Jet_pt_11_ = Jet_pt_12_ = Jet_pt_13_ = Jet_pt_14_ = Jet_pt_15_ = Jet_pt_16_ = Jet_pt_17_ = Jet_pt_18_ = Jet_pt_19_ = Jet_pt_20_ = -1;
        Jet_eta_1_ = Jet_eta_2_ = Jet_eta_3_ = Jet_eta_4_ = Jet_eta_5_ = Jet_eta_6_ = Jet_eta_7_ = Jet_eta_8_ = Jet_eta_9_ = Jet_eta_10_ = -1;
        Jet_eta_11_ = Jet_eta_12_ = Jet_eta_13_ = Jet_eta_14_ = Jet_eta_15_ = Jet_eta_16_ = Jet_eta_17_ = Jet_eta_18_ = Jet_eta_19_ = Jet_eta_20_ = -1;
        Jet_phi_1_ = Jet_phi_2_ = Jet_phi_3_ = Jet_phi_4_ = Jet_phi_5_ = Jet_phi_6_ = Jet_phi_7_ = Jet_phi_8_ = Jet_phi_9_ = Jet_phi_10_ = -1;
        Jet_phi_11_ = Jet_phi_12_ = Jet_phi_13_ = Jet_phi_14_ = Jet_phi_15_ = Jet_phi_16_ = Jet_phi_17_ = Jet_phi_18_ = Jet_phi_19_ = Jet_phi_20_ = -1;
        Jet_m_1_ = Jet_m_2_ = Jet_m_3_ = Jet_m_4_ = Jet_m_5_ = Jet_m_6_ = Jet_m_7_ = Jet_m_8_ = Jet_m_9_ = Jet_m_10_ = -1;
        Jet_m_11_ = Jet_m_12_ = Jet_m_13_ = Jet_m_14_ = Jet_m_15_ = Jet_m_16_ = Jet_m_17_ = Jet_m_18_ = Jet_m_19_ = Jet_m_20_ = -1;
        GoodLeptons_pt_1_ = GoodLeptons_eta_1_ = GoodLeptons_phi_1_ = GoodLeptons_m_1_ = -1;
        BestComboAvgMass_ = -1;
    }

    /**
     *The job of mapVars is to populate the internal offests for all variables in the input variable list with their memory location in the data array.  To be called only once.
     */
    void mapVars(const std::vector<std::string>& vars)
    {
        len_ = vars.size();

        for(unsigned int i = 0; i < vars.size(); ++i)
        {
            if(     vars[i].compare("fwm2_top6") == 0)  fwm2_top6_ = i;
            else if(vars[i].compare("fwm3_top6") == 0)  fwm3_top6_ = i;
            else if(vars[i].compare("fwm4_top6") == 0)  fwm4_top6_ = i;
            else if(vars[i].compare("fwm5_top6") == 0)  fwm5_top6_ = i;
            else if(vars[i].compare("fwm6_top6") == 0)  fwm6_top6_ = i;
            else if(vars[i].compare("fwm7_top6") == 0)  fwm7_top6_ = i;
            else if(vars[i].compare("fwm8_top6") == 0)  fwm8_top6_ = i;
            else if(vars[i].compare("fwm9_top6") == 0)  fwm9_top6_ = i;
            else if(vars[i].compare("fwm10_top6") == 0) fwm10_top6_ = i;
            else if(vars[i].compare("jmt_ev0_top6") == 0) jmt_ev0_top6_ = i;
            else if(vars[i].compare("jmt_ev1_top6") == 0) jmt_ev1_top6_ = i;
            else if(vars[i].compare("jmt_ev2_top6") == 0) jmt_ev2_top6_ = i;
            else if(vars[i].compare("NGoodJets_double") == 0) NGoodJets_double_ = i;
            else if(vars[i].compare("Jet_pt_1") == 0) Jet_pt_1_ = i;
            else if(vars[i].compare("Jet_pt_2") == 0) Jet_pt_2_ = i;
            else if(vars[i].compare("Jet_pt_3") == 0) Jet_pt_3_ = i;
            else if(vars[i].compare("Jet_pt_4") == 0) Jet_pt_4_ = i;
            else if(vars[i].compare("Jet_pt_5") == 0) Jet_pt_5_ = i;
            else if(vars[i].compare("Jet_pt_6") == 0) Jet_pt_6_ = i;
            else if(vars[i].compare("Jet_pt_7") == 0) Jet_pt_7_ = i;
            else if(vars[i].compare("Jet_pt_8") == 0) Jet_pt_8_ = i;
            else if(vars[i].compare("Jet_pt_9") == 0) Jet_pt_9_ = i;
            else if(vars[i].compare("Jet_pt_10") == 0) Jet_pt_10_ = i;
            else if(vars[i].compare("Jet_pt_11") == 0) Jet_pt_11_ = i;
            else if(vars[i].compare("Jet_pt_12") == 0) Jet_pt_12_ = i;
            else if(vars[i].compare("Jet_pt_13") == 0) Jet_pt_13_ = i;
            else if(vars[i].compare("Jet_pt_14") == 0) Jet_pt_14_ = i;
            else if(vars[i].compare("Jet_pt_15") == 0) Jet_pt_15_ = i;
            else if(vars[i].compare("Jet_pt_16") == 0) Jet_pt_16_ = i;
            else if(vars[i].compare("Jet_pt_17") == 0) Jet_pt_17_ = i;
            else if(vars[i].compare("Jet_pt_18") == 0) Jet_pt_18_ = i;
            else if(vars[i].compare("Jet_pt_19") == 0) Jet_pt_19_ = i;
            else if(vars[i].compare("Jet_pt_20") == 0) Jet_pt_20_ = i;
            else if(vars[i].compare("Jet_eta_1") == 0) Jet_eta_1_ = i;
            else if(vars[i].compare("Jet_eta_2") == 0) Jet_eta_2_ = i;
            else if(vars[i].compare("Jet_eta_3") == 0) Jet_eta_3_ = i;
            else if(vars[i].compare("Jet_eta_4") == 0) Jet_eta_4_ = i;
            else if(vars[i].compare("Jet_eta_5") == 0) Jet_eta_5_ = i;
            else if(vars[i].compare("Jet_eta_6") == 0) Jet_eta_6_ = i;
            else if(vars[i].compare("Jet_eta_7") == 0) Jet_eta_7_ = i;
            else if(vars[i].compare("Jet_eta_8") == 0) Jet_eta_8_ = i;
            else if(vars[i].compare("Jet_eta_9") == 0) Jet_eta_9_ = i;
            else if(vars[i].compare("Jet_eta_10") == 0) Jet_eta_10_ = i;
            else if(vars[i].compare("Jet_eta_11") == 0) Jet_eta_11_ = i;
            else if(vars[i].compare("Jet_eta_12") == 0) Jet_eta_12_ = i;
            else if(vars[i].compare("Jet_eta_13") == 0) Jet_eta_13_ = i;
            else if(vars[i].compare("Jet_eta_14") == 0) Jet_eta_14_ = i;
            else if(vars[i].compare("Jet_eta_15") == 0) Jet_eta_15_ = i;
            else if(vars[i].compare("Jet_eta_16") == 0) Jet_eta_16_ = i;
            else if(vars[i].compare("Jet_eta_17") == 0) Jet_eta_17_ = i;
            else if(vars[i].compare("Jet_eta_18") == 0) Jet_eta_18_ = i;
            else if(vars[i].compare("Jet_eta_19") == 0) Jet_eta_19_ = i;
            else if(vars[i].compare("Jet_eta_20") == 0) Jet_eta_20_ = i;
            else if(vars[i].compare("Jet_phi_1") == 0) Jet_phi_1_ = i;
            else if(vars[i].compare("Jet_phi_2") == 0) Jet_phi_2_ = i;
            else if(vars[i].compare("Jet_phi_3") == 0) Jet_phi_3_ = i;
            else if(vars[i].compare("Jet_phi_4") == 0) Jet_phi_4_ = i;
            else if(vars[i].compare("Jet_phi_5") == 0) Jet_phi_5_ = i;
            else if(vars[i].compare("Jet_phi_6") == 0) Jet_phi_6_ = i;
            else if(vars[i].compare("Jet_phi_7") == 0) Jet_phi_7_ = i;
            else if(vars[i].compare("Jet_phi_8") == 0) Jet_phi_8_ = i;
            else if(vars[i].compare("Jet_phi_9") == 0) Jet_phi_9_ = i;
            else if(vars[i].compare("Jet_phi_10") == 0) Jet_phi_10_ = i;
            else if(vars[i].compare("Jet_phi_11") == 0) Jet_phi_11_ = i;
            else if(vars[i].compare("Jet_phi_12") == 0) Jet_phi_12_ = i;
            else if(vars[i].compare("Jet_phi_13") == 0) Jet_phi_13_ = i;
            else if(vars[i].compare("Jet_phi_14") == 0) Jet_phi_14_ = i;
            else if(vars[i].compare("Jet_phi_15") == 0) Jet_phi_15_ = i;
            else if(vars[i].compare("Jet_phi_16") == 0) Jet_phi_16_ = i;
            else if(vars[i].compare("Jet_phi_17") == 0) Jet_phi_17_ = i;
            else if(vars[i].compare("Jet_phi_18") == 0) Jet_phi_18_ = i;
            else if(vars[i].compare("Jet_phi_19") == 0) Jet_phi_19_ = i;
            else if(vars[i].compare("Jet_phi_20") == 0) Jet_phi_20_ = i;
            else if(vars[i].compare("Jet_m_1") == 0) Jet_m_1_ = i;
            else if(vars[i].compare("Jet_m_2") == 0) Jet_m_2_ = i;
            else if(vars[i].compare("Jet_m_3") == 0) Jet_m_3_ = i;
            else if(vars[i].compare("Jet_m_4") == 0) Jet_m_4_ = i;
            else if(vars[i].compare("Jet_m_5") == 0) Jet_m_5_ = i;
            else if(vars[i].compare("Jet_m_6") == 0) Jet_m_6_ = i;
            else if(vars[i].compare("Jet_m_7") == 0) Jet_m_7_ = i;
            else if(vars[i].compare("Jet_m_8") == 0) Jet_m_8_ = i;
            else if(vars[i].compare("Jet_m_9") == 0) Jet_m_9_ = i;
            else if(vars[i].compare("Jet_m_10") == 0) Jet_m_10_ = i;
            else if(vars[i].compare("Jet_m_11") == 0) Jet_m_11_ = i;
            else if(vars[i].compare("Jet_m_12") == 0) Jet_m_12_ = i;
            else if(vars[i].compare("Jet_m_13") == 0) Jet_m_13_ = i;
            else if(vars[i].compare("Jet_m_14") == 0) Jet_m_14_ = i;
            else if(vars[i].compare("Jet_m_15") == 0) Jet_m_15_ = i;
            else if(vars[i].compare("Jet_m_16") == 0) Jet_m_16_ = i;
            else if(vars[i].compare("Jet_m_17") == 0) Jet_m_17_ = i;
            else if(vars[i].compare("Jet_m_18") == 0) Jet_m_18_ = i;
            else if(vars[i].compare("Jet_m_19") == 0) Jet_m_19_ = i;
            else if(vars[i].compare("Jet_m_20") == 0) Jet_m_20_ = i;
            else if(vars[i].compare("GoodLeptons_pt_1") == 0) GoodLeptons_pt_1_ = i;
            else if(vars[i].compare("GoodLeptons_eta_1") == 0) GoodLeptons_eta_1_ = i;
            else if(vars[i].compare("GoodLeptons_phi_1") == 0) GoodLeptons_phi_1_ = i;
            else if(vars[i].compare("GoodLeptons_m_1") == 0) GoodLeptons_m_1_ = i;
            else if(vars[i].compare("BestComboAvgMass") == 0) BestComboAvgMass_ = i;
        }
    }
    /**
     *The job of setPtr is to set the starting place of memory block where the data will be written. To be called only once for the creation of the array pointed to by data.
     */
    void setPtr(float* data) {basePtr_ = data;}
    /**
     *Calculate the requested variables and store the values directly in the input array for the MVA
     */
    void calculateVars(const NTupleReader& tr, int iCand)
    {
        if(fwm2_top6_ >= 0)  *(basePtr_ + fwm2_top6_ + len_*iCand) =  tr.getVar<double>("fwm2_top6");
        if(fwm3_top6_ >= 0)  *(basePtr_ + fwm3_top6_ + len_*iCand) =  tr.getVar<double>("fwm3_top6");
        if(fwm4_top6_ >= 0)  *(basePtr_ + fwm4_top6_ + len_*iCand) =  tr.getVar<double>("fwm4_top6");
        if(fwm5_top6_ >= 0)  *(basePtr_ + fwm5_top6_ + len_*iCand) =  tr.getVar<double>("fwm5_top6");
        if(fwm6_top6_ >= 0)  *(basePtr_ + fwm6_top6_ + len_*iCand) =  tr.getVar<double>("fwm6_top6");
        if(fwm7_top6_ >= 0)  *(basePtr_ + fwm7_top6_ + len_*iCand) =  tr.getVar<double>("fwm7_top6");
        if(fwm8_top6_ >= 0)  *(basePtr_ + fwm8_top6_ + len_*iCand) =  tr.getVar<double>("fwm8_top6");
        if(fwm9_top6_ >= 0)  *(basePtr_ + fwm9_top6_ + len_*iCand) =  tr.getVar<double>("fwm9_top6");
        if(fwm10_top6_ >= 0) *(basePtr_ + fwm10_top6_ + len_*iCand) =  tr.getVar<double>("fwm10_top6");
        if(jmt_ev0_top6_ >= 0) *(basePtr_ + jmt_ev0_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev0_top6");
        if(jmt_ev1_top6_ >= 0) *(basePtr_ + jmt_ev1_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev1_top6");
        if(jmt_ev2_top6_ >= 0) *(basePtr_ + jmt_ev2_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev2_top6");
        if(NGoodJets_double_ >= 0) *(basePtr_ + NGoodJets_double_ + len_*iCand) =  static_cast<double>(tr.getVar<unsigned long>("NGoodJets"));
        if(Jet_pt_1_  >= 0) *(basePtr_ + Jet_pt_1_  + len_*iCand) = tr.getVar<double>("Jet_pt_1");
        if(Jet_pt_2_  >= 0) *(basePtr_ + Jet_pt_2_  + len_*iCand) = tr.getVar<double>("Jet_pt_2");
        if(Jet_pt_3_  >= 0) *(basePtr_ + Jet_pt_3_  + len_*iCand) = tr.getVar<double>("Jet_pt_3");
        if(Jet_pt_4_  >= 0) *(basePtr_ + Jet_pt_4_  + len_*iCand) = tr.getVar<double>("Jet_pt_4");
        if(Jet_pt_5_  >= 0) *(basePtr_ + Jet_pt_5_  + len_*iCand) = tr.getVar<double>("Jet_pt_5");
        if(Jet_pt_6_  >= 0) *(basePtr_ + Jet_pt_6_  + len_*iCand) = tr.getVar<double>("Jet_pt_6");
        if(Jet_pt_7_  >= 0) *(basePtr_ + Jet_pt_7_  + len_*iCand) = tr.getVar<double>("Jet_pt_7");
        if(Jet_pt_8_  >= 0) *(basePtr_ + Jet_pt_8_  + len_*iCand) = tr.getVar<double>("Jet_pt_8");
        if(Jet_pt_9_  >= 0) *(basePtr_ + Jet_pt_9_  + len_*iCand) = tr.getVar<double>("Jet_pt_9");
        if(Jet_pt_10_ >= 0) *(basePtr_ + Jet_pt_10_ + len_*iCand) = tr.getVar<double>("Jet_pt_10");
        if(Jet_pt_11_ >= 0) *(basePtr_ + Jet_pt_11_ + len_*iCand) = tr.getVar<double>("Jet_pt_11");
        if(Jet_pt_12_ >= 0) *(basePtr_ + Jet_pt_12_ + len_*iCand) = tr.getVar<double>("Jet_pt_12");
        if(Jet_pt_13_ >= 0) *(basePtr_ + Jet_pt_13_ + len_*iCand) = tr.getVar<double>("Jet_pt_13");
        if(Jet_pt_14_ >= 0) *(basePtr_ + Jet_pt_14_ + len_*iCand) = tr.getVar<double>("Jet_pt_14");
        if(Jet_pt_15_ >= 0) *(basePtr_ + Jet_pt_15_ + len_*iCand) = tr.getVar<double>("Jet_pt_15");
        if(Jet_pt_16_ >= 0) *(basePtr_ + Jet_pt_16_ + len_*iCand) = tr.getVar<double>("Jet_pt_16");
        if(Jet_pt_17_ >= 0) *(basePtr_ + Jet_pt_17_ + len_*iCand) = tr.getVar<double>("Jet_pt_17");
        if(Jet_pt_18_ >= 0) *(basePtr_ + Jet_pt_18_ + len_*iCand) = tr.getVar<double>("Jet_pt_18");
        if(Jet_pt_19_ >= 0) *(basePtr_ + Jet_pt_19_ + len_*iCand) = tr.getVar<double>("Jet_pt_19");
        if(Jet_pt_20_ >= 0) *(basePtr_ + Jet_pt_20_ + len_*iCand) = tr.getVar<double>("Jet_pt_20");
        if(Jet_eta_1_  >= 0) *(basePtr_ + Jet_eta_1_  + len_*iCand) = tr.getVar<double>("Jet_eta_1");
        if(Jet_eta_2_  >= 0) *(basePtr_ + Jet_eta_2_  + len_*iCand) = tr.getVar<double>("Jet_eta_2");
        if(Jet_eta_3_  >= 0) *(basePtr_ + Jet_eta_3_  + len_*iCand) = tr.getVar<double>("Jet_eta_3");
        if(Jet_eta_4_  >= 0) *(basePtr_ + Jet_eta_4_  + len_*iCand) = tr.getVar<double>("Jet_eta_4");
        if(Jet_eta_5_  >= 0) *(basePtr_ + Jet_eta_5_  + len_*iCand) = tr.getVar<double>("Jet_eta_5");
        if(Jet_eta_6_  >= 0) *(basePtr_ + Jet_eta_6_  + len_*iCand) = tr.getVar<double>("Jet_eta_6");
        if(Jet_eta_7_  >= 0) *(basePtr_ + Jet_eta_7_  + len_*iCand) = tr.getVar<double>("Jet_eta_7");
        if(Jet_eta_8_  >= 0) *(basePtr_ + Jet_eta_8_  + len_*iCand) = tr.getVar<double>("Jet_eta_8");
        if(Jet_eta_9_  >= 0) *(basePtr_ + Jet_eta_9_  + len_*iCand) = tr.getVar<double>("Jet_eta_9");
        if(Jet_eta_10_ >= 0) *(basePtr_ + Jet_eta_10_ + len_*iCand) = tr.getVar<double>("Jet_eta_10");
        if(Jet_eta_11_ >= 0) *(basePtr_ + Jet_eta_11_ + len_*iCand) = tr.getVar<double>("Jet_eta_11");
        if(Jet_eta_12_ >= 0) *(basePtr_ + Jet_eta_12_ + len_*iCand) = tr.getVar<double>("Jet_eta_12");
        if(Jet_eta_13_ >= 0) *(basePtr_ + Jet_eta_13_ + len_*iCand) = tr.getVar<double>("Jet_eta_13");
        if(Jet_eta_14_ >= 0) *(basePtr_ + Jet_eta_14_ + len_*iCand) = tr.getVar<double>("Jet_eta_14");
        if(Jet_eta_15_ >= 0) *(basePtr_ + Jet_eta_15_ + len_*iCand) = tr.getVar<double>("Jet_eta_15");
        if(Jet_eta_16_ >= 0) *(basePtr_ + Jet_eta_16_ + len_*iCand) = tr.getVar<double>("Jet_eta_16");
        if(Jet_eta_17_ >= 0) *(basePtr_ + Jet_eta_17_ + len_*iCand) = tr.getVar<double>("Jet_eta_17");
        if(Jet_eta_18_ >= 0) *(basePtr_ + Jet_eta_18_ + len_*iCand) = tr.getVar<double>("Jet_eta_18");
        if(Jet_eta_19_ >= 0) *(basePtr_ + Jet_eta_19_ + len_*iCand) = tr.getVar<double>("Jet_eta_19");
        if(Jet_eta_20_ >= 0) *(basePtr_ + Jet_eta_20_ + len_*iCand) = tr.getVar<double>("Jet_eta_20");
        if(Jet_phi_1_  >= 0) *(basePtr_ + Jet_phi_1_  + len_*iCand) = tr.getVar<double>("Jet_phi_1");
        if(Jet_phi_2_  >= 0) *(basePtr_ + Jet_phi_2_  + len_*iCand) = tr.getVar<double>("Jet_phi_2");
        if(Jet_phi_3_  >= 0) *(basePtr_ + Jet_phi_3_  + len_*iCand) = tr.getVar<double>("Jet_phi_3");
        if(Jet_phi_4_  >= 0) *(basePtr_ + Jet_phi_4_  + len_*iCand) = tr.getVar<double>("Jet_phi_4");
        if(Jet_phi_5_  >= 0) *(basePtr_ + Jet_phi_5_  + len_*iCand) = tr.getVar<double>("Jet_phi_5");
        if(Jet_phi_6_  >= 0) *(basePtr_ + Jet_phi_6_  + len_*iCand) = tr.getVar<double>("Jet_phi_6");
        if(Jet_phi_7_  >= 0) *(basePtr_ + Jet_phi_7_  + len_*iCand) = tr.getVar<double>("Jet_phi_7");
        if(Jet_phi_8_  >= 0) *(basePtr_ + Jet_phi_8_  + len_*iCand) = tr.getVar<double>("Jet_phi_8");
        if(Jet_phi_9_  >= 0) *(basePtr_ + Jet_phi_9_  + len_*iCand) = tr.getVar<double>("Jet_phi_9");
        if(Jet_phi_10_ >= 0) *(basePtr_ + Jet_phi_10_ + len_*iCand) = tr.getVar<double>("Jet_phi_10");
        if(Jet_phi_11_ >= 0) *(basePtr_ + Jet_phi_11_ + len_*iCand) = tr.getVar<double>("Jet_phi_11");
        if(Jet_phi_12_ >= 0) *(basePtr_ + Jet_phi_12_ + len_*iCand) = tr.getVar<double>("Jet_phi_12");
        if(Jet_phi_13_ >= 0) *(basePtr_ + Jet_phi_13_ + len_*iCand) = tr.getVar<double>("Jet_phi_13");
        if(Jet_phi_14_ >= 0) *(basePtr_ + Jet_phi_14_ + len_*iCand) = tr.getVar<double>("Jet_phi_14");
        if(Jet_phi_15_ >= 0) *(basePtr_ + Jet_phi_15_ + len_*iCand) = tr.getVar<double>("Jet_phi_15");
        if(Jet_phi_16_ >= 0) *(basePtr_ + Jet_phi_16_ + len_*iCand) = tr.getVar<double>("Jet_phi_16");
        if(Jet_phi_17_ >= 0) *(basePtr_ + Jet_phi_17_ + len_*iCand) = tr.getVar<double>("Jet_phi_17");
        if(Jet_phi_18_ >= 0) *(basePtr_ + Jet_phi_18_ + len_*iCand) = tr.getVar<double>("Jet_phi_18");
        if(Jet_phi_19_ >= 0) *(basePtr_ + Jet_phi_19_ + len_*iCand) = tr.getVar<double>("Jet_phi_19");
        if(Jet_phi_20_ >= 0) *(basePtr_ + Jet_phi_20_ + len_*iCand) = tr.getVar<double>("Jet_phi_20");
        if(Jet_m_1_  >= 0) *(basePtr_ + Jet_m_1_  + len_*iCand) = tr.getVar<double>("Jet_m_1");
        if(Jet_m_2_  >= 0) *(basePtr_ + Jet_m_2_  + len_*iCand) = tr.getVar<double>("Jet_m_2");
        if(Jet_m_3_  >= 0) *(basePtr_ + Jet_m_3_  + len_*iCand) = tr.getVar<double>("Jet_m_3");
        if(Jet_m_4_  >= 0) *(basePtr_ + Jet_m_4_  + len_*iCand) = tr.getVar<double>("Jet_m_4");
        if(Jet_m_5_  >= 0) *(basePtr_ + Jet_m_5_  + len_*iCand) = tr.getVar<double>("Jet_m_5");
        if(Jet_m_6_  >= 0) *(basePtr_ + Jet_m_6_  + len_*iCand) = tr.getVar<double>("Jet_m_6");
        if(Jet_m_7_  >= 0) *(basePtr_ + Jet_m_7_  + len_*iCand) = tr.getVar<double>("Jet_m_7");
        if(Jet_m_8_  >= 0) *(basePtr_ + Jet_m_8_  + len_*iCand) = tr.getVar<double>("Jet_m_8");
        if(Jet_m_9_  >= 0) *(basePtr_ + Jet_m_9_  + len_*iCand) = tr.getVar<double>("Jet_m_9");
        if(Jet_m_10_ >= 0) *(basePtr_ + Jet_m_10_ + len_*iCand) = tr.getVar<double>("Jet_m_10");
        if(Jet_m_11_ >= 0) *(basePtr_ + Jet_m_11_ + len_*iCand) = tr.getVar<double>("Jet_m_11");
        if(Jet_m_12_ >= 0) *(basePtr_ + Jet_m_12_ + len_*iCand) = tr.getVar<double>("Jet_m_12");
        if(Jet_m_13_ >= 0) *(basePtr_ + Jet_m_13_ + len_*iCand) = tr.getVar<double>("Jet_m_13");
        if(Jet_m_14_ >= 0) *(basePtr_ + Jet_m_14_ + len_*iCand) = tr.getVar<double>("Jet_m_14");
        if(Jet_m_15_ >= 0) *(basePtr_ + Jet_m_15_ + len_*iCand) = tr.getVar<double>("Jet_m_15");
        if(Jet_m_16_ >= 0) *(basePtr_ + Jet_m_16_ + len_*iCand) = tr.getVar<double>("Jet_m_16");
        if(Jet_m_17_ >= 0) *(basePtr_ + Jet_m_17_ + len_*iCand) = tr.getVar<double>("Jet_m_17");
        if(Jet_m_18_ >= 0) *(basePtr_ + Jet_m_18_ + len_*iCand) = tr.getVar<double>("Jet_m_18");
        if(Jet_m_19_ >= 0) *(basePtr_ + Jet_m_19_ + len_*iCand) = tr.getVar<double>("Jet_m_19");
        if(Jet_m_20_ >= 0) *(basePtr_ + Jet_m_20_ + len_*iCand) = tr.getVar<double>("Jet_m_20");
        if(GoodLeptons_pt_1_  >= 0) *(basePtr_ + GoodLeptons_pt_1_  + len_*iCand) = tr.getVar<double>("GoodLeptons_pt_1");
        if(GoodLeptons_eta_1_ >= 0) *(basePtr_ + GoodLeptons_eta_1_ + len_*iCand) = tr.getVar<double>("GoodLeptons_eta_1");
        if(GoodLeptons_phi_1_ >= 0) *(basePtr_ + GoodLeptons_phi_1_ + len_*iCand) = tr.getVar<double>("GoodLeptons_phi_1");
        if(GoodLeptons_m_1_   >= 0) *(basePtr_ + GoodLeptons_m_1_   + len_*iCand) = tr.getVar<double>("GoodLeptons_m_1");
        if(BestComboAvgMass_ >= 0) *(basePtr_ + BestComboAvgMass_ + len_*iCand) = tr.getVar<double>("BestComboAvgMass");
    }
};

class DeepEventShape
{
private:
    double discriminator_;
    std::string modelFile_, inputOp_, outputOp_, nJetMask_;

    //Tensoflow session pointer
    TF_Session* session_;

    //Input variable names 
    std::vector<std::string> vars_;
    std::vector<double> binEdges_;

    std::vector<TF_Output>     inputs_;
    std::vector<TF_Output>     outputs_;
    std::vector<TF_Operation*> targets_;

    //variable calclator
    std::shared_ptr<EventShapeCalculator> varCalculator_;

    template<typename T> std::vector<T> getVecFromCfg(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& var, const cfg::Context& localCxt, const T& defaultYo)
    {
        std::vector<T> vec;
        int iVar = 0;
        bool keepLooping;
        do
        {
            keepLooping = false;
        
            //Get variable name
            T v = cfgDoc->get(var, iVar, localCxt, defaultYo);
        
            //if it is a non empty string save in vector
            if(v != defaultYo)
            {
                keepLooping = true;
        
                vec.push_back(v);
            }
            ++iVar;
        }
        while(keepLooping);
        
        return vec;
    }

    void getParameters(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& localContextName)
    {
        //Construct contexts
        cfg::Context localCxt(localContextName);

        modelFile_     = cfgDoc->get("modelFile", localCxt, "");
        inputOp_       = cfgDoc->get("inputOp",   localCxt, "x");
        outputOp_      = cfgDoc->get("outputOp",  localCxt, "y");
        nJetMask_      = cfgDoc->get("nJetMask",  localCxt, "");
        vars_          = getVecFromCfg<std::string>(cfgDoc, "mvaVar", localCxt, "");
        binEdges_      = getVecFromCfg<double>(cfgDoc, "binEdges", localCxt, -1);
        
        //Variable to hold tensorflow status
        TF_Status* status = TF_NewStatus();
        
        //get the grafdef from the file
        TF_Buffer* graph_def = read_file(modelFile_);
        
        // Import graph_def into graph
        TF_Graph* graph = TF_NewGraph();
        TF_ImportGraphDefOptions* graph_opts = TF_NewImportGraphDefOptions();
        TF_GraphImportGraphDef(graph, graph_def, graph_opts, status);
        TF_DeleteImportGraphDefOptions(graph_opts);
        TF_DeleteBuffer(graph_def);
        
        //Create tensorflow session from imported graph
        TF_SessionOptions* sess_opts = TF_NewSessionOptions();
        uint8_t config[] = {0x10, 0x01};
        TF_SetConfig(sess_opts, static_cast<void*>(config), 2, status);
        session_ = TF_NewSession(graph, sess_opts, status);
        TF_DeleteSessionOptions(sess_opts);
        
        TF_Operation* op_x = TF_GraphOperationByName(graph, inputOp_.c_str());
        TF_Operation* op_y = TF_GraphOperationByName(graph, outputOp_.c_str());
        
        //Clean up graph
        TF_DeleteGraph(graph);
        
        inputs_ .emplace_back(TF_Output({op_x, 0}));
        outputs_.emplace_back(TF_Output({op_y, 0}));
        targets_.emplace_back(op_y);
        
        TF_DeleteStatus(status);

        //map variables
        varCalculator_.reset(new EventShapeCalculator());
        varCalculator_->mapVars(vars_);
    }

    void runDeepEventShape(NTupleReader& tr)
    {
        //tensorflow status variable
        TF_Status* status = TF_NewStatus();
        
        //Create place to store the output vectors 
        std::vector<TF_Tensor*>    output_values(1);
        
        //Construct tensorflow input tensor
        std::vector<TF_Tensor*> input_values;
        const int elemSize = sizeof(float);
        std::vector<int64_t> dims = {static_cast<int64_t>(1), static_cast<int64_t>(vars_.size())};
        int nelem = 1;
        for(const auto dimLen : dims) nelem *= dimLen;
        TF_Tensor* input_values_0 =  TF_AllocateTensor(TF_FLOAT, dims.data(), dims.size(), elemSize*nelem);
        
        input_values = { input_values_0 };
        varCalculator_->setPtr(static_cast<float*>(TF_TensorData(input_values_0)));

        int iCand = 0;
        varCalculator_->calculateVars(tr, iCand);

        //predict values
        TF_SessionRun(session_,
                      // RunOptions
                      nullptr,
                      // Input tensors
                      inputs_.data(), input_values.data(), inputs_.size(),
                      // Output tensors
                      outputs_.data(), output_values.data(), outputs_.size(),
                      // Target operations
                      targets_.data(), targets_.size(),
                      // RunMetadata
                      nullptr,
                      // Output status
                      status);
        
        //Get output discriminators 
        auto discriminators = static_cast<float*>(TF_TensorData(output_values[0]));                
        
        //discriminators is a 2D array, we only want the first entry of every array
        double discriminator = static_cast<double>(discriminators[iCand*TF_Dim(output_values[0], 1)]);

        for(auto tensor : input_values)  TF_DeleteTensor(tensor);
        for(auto tensor : output_values) TF_DeleteTensor(tensor);
        
        TF_DeleteStatus(status);

        // Register Variables
        tr.registerDerivedVar("deepESM"+nJetMask_+"_val", discriminator);
        // Define and register deepESM bins
        for(int i = 1; i < binEdges_.size(); i++)
        {
            bool passDeepESMBin = discriminator > binEdges_[i-1] && discriminator <= binEdges_[i];
            tr.registerDerivedVar("deepESM"+nJetMask_+"_bin"+std::to_string(i), passDeepESMBin);
        }
    }

    static void free_buffer(void* data, size_t length) 
    {
        free(data);
    }

    TF_Buffer* read_file(const std::string& file) 
    {
        FILE* f = fopen(file.c_str(), "rb");

        fseek(f, 0, SEEK_END);
        long fsize = ftell(f);
        fseek(f, 0, SEEK_SET);  //same as rewind(f);

        void* data = malloc(fsize);
        fread(data, fsize, 1, f);
        fclose(f);

        TF_Buffer* buf = TF_NewBuffer();
        buf->data = data;
        buf->length = fsize;
        buf->data_deallocator = free_buffer;
        return buf;
    }

public:
    DeepEventShape(const std::string cfgFileName = "DeepEventShape.cfg", std::string localContextName = "Info")
    {
        std::cout<<"Setting up DeepEventShape"<<std::endl;
        //buffer to hold file contents 
        std::string cfgText;

        FILE *f = fopen(cfgFileName.c_str(), "r");
        char buff[1024];
        for(; !feof(f) && fgets(buff, 1023, f);)
        {
            cfgText += buff;
        }
        
        fclose(f);
        
        //pass raw text to cfg parser, to return parsed document
        std::unique_ptr<cfg::CfgDocument> cfgDoc = cfg::CfgDocument::parseDocument(cfgText);
        getParameters(cfgDoc, localContextName);

        std::cout<<"Using "+cfgFileName+" as the DeepEventShape config file"<<std::endl;
    }

    ~DeepEventShape()
    {
        //tensorflow status variable
        //TF_Status* status = TF_NewStatus();
        //TF_DeleteSession(session_, status);
        //TF_DeleteStatus(status);
    }
    
    void operator()(NTupleReader& tr)
    {
        runDeepEventShape(tr);
    }
};

#endif
