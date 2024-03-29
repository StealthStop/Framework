#ifndef PREPNTUPLEVARS_H
#define PREPNTUPLEVARS_H

#include "Framework/Framework/include/Utility.h"

class PrepNTupleVars
{
private:
    std::map<std::string, std::pair<double,double>> pTMass_;

    class JetCollection
    {
    public:
        const std::vector<utility::LorentzVector>& Jets;
        const std::vector<TLorentzVector>&         Jets_TLV;
        const std::vector<float>&                  Jets_bDiscriminatorCSV;
        const std::vector<float>&                  Jets_bJetTagDeepCSVprobb;
        const std::vector<float>&                  Jets_bJetTagDeepCSVprobbb;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourprobb;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourprobbb;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourproblepb;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourprobc;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourprobuds;
        const std::vector<float>&                  Jets_bJetTagDeepFlavourprobg;
        const std::vector<float>&                  Jets_bJetTagDeepCSVprobc;
        const std::vector<float>&                  Jets_bJetTagDeepCSVprobudsg;
        const std::vector<float>&                  Jets_qgLikelihood;
        const std::vector<float>&                  Jets_muonEnergyFraction;
        const std::vector<float>&                  Jets_hfHadronEnergyFraction;
        const std::vector<float>&                  Jets_hfEMEnergyFraction;
        const std::vector<float>&                  Jets_photonEnergyFraction;
        const std::vector<float>&                  Jets_electronEnergyFraction;
        const std::vector<int>&                    Jets_chargedHadronMultiplicity;
        const std::vector<int>&                    Jets_neutralHadronMultiplicity;
        const std::vector<int>&                    Jets_photonMultiplicity;
        const std::vector<int>&                    Jets_electronMultiplicity;
        const std::vector<int>&                    Jets_muonMultiplicity;
        const std::vector<bool>&                   Jets_ID;
        const bool&                                JetID;
        const std::vector<int>&                    Jets_partonFlavor;
        const std::vector<float>&                  Jets_ptD;
        const std::vector<float>&                  Jets_axismajor;
        const std::vector<float>&                  Jets_axisminor;
        const std::vector<int>&                    Jets_multiplicity;
        const std::vector<float>&                  Jets_neutralEmEnergyFraction;
        const std::vector<float>&                  Jets_chargedEmEnergyFraction;
        const std::vector<float>&                  Jets_neutralHadronEnergyFraction;
        const std::vector<float>&                  Jets_chargedHadronEnergyFraction;
        const std::vector<float>&                  Jets_pileupId;
        
    JetCollection(const NTupleReader& tr) 
            : Jets(tr.getVec<utility::LorentzVector>("Jets"))
            , Jets_TLV(tr.getVec<TLorentzVector>("Jets_TLV"))
            , Jets_bDiscriminatorCSV(tr.getVec<float>("Jets_bDiscriminatorCSV")) 
            , Jets_bJetTagDeepCSVprobb(tr.getVec<float>("Jets_bJetTagDeepCSVprobb"))
            , Jets_bJetTagDeepCSVprobbb(tr.getVec<float>("Jets_bJetTagDeepCSVprobbb"))
            , Jets_bJetTagDeepFlavourprobb(tr.getVec<float>("Jets_bJetTagDeepFlavourprobb"))
            , Jets_bJetTagDeepFlavourprobbb(tr.getVec<float>("Jets_bJetTagDeepFlavourprobbb"))
            , Jets_bJetTagDeepFlavourproblepb(tr.getVec<float>("Jets_bJetTagDeepFlavourproblepb"))
            , Jets_bJetTagDeepFlavourprobc(tr.getVec<float>("Jets_bJetTagDeepFlavourprobc"))
            , Jets_bJetTagDeepFlavourprobuds(tr.getVec<float>("Jets_bJetTagDeepFlavourprobuds"))
            , Jets_bJetTagDeepFlavourprobg(tr.getVec<float>("Jets_bJetTagDeepFlavourprobg"))
            , Jets_bJetTagDeepCSVprobc(tr.getVec<float>("Jets_bJetTagDeepCSVprobc"))
            , Jets_bJetTagDeepCSVprobudsg(tr.getVec<float>("Jets_bJetTagDeepCSVprobudsg"))
            , Jets_qgLikelihood(tr.getVec<float>("Jets_qgLikelihood"))
            , Jets_muonEnergyFraction(tr.getVec<float>("Jets_muonEnergyFraction"))
            , Jets_hfHadronEnergyFraction(tr.getVec<float>("Jets_hfHadronEnergyFraction"))
            , Jets_hfEMEnergyFraction(tr.getVec<float>("Jets_hfEMEnergyFraction"))
            , Jets_photonEnergyFraction(tr.getVec<float>("Jets_photonEnergyFraction"))
            , Jets_electronEnergyFraction(tr.getVec<float>("Jets_electronEnergyFraction"))
            , Jets_chargedHadronMultiplicity(tr.getVec<int>("Jets_chargedHadronMultiplicity"))
            , Jets_neutralHadronMultiplicity(tr.getVec<int>("Jets_neutralHadronMultiplicity"))
            , Jets_photonMultiplicity(tr.getVec<int>("Jets_photonMultiplicity"))
            , Jets_electronMultiplicity(tr.getVec<int>("Jets_electronMultiplicity"))
            , Jets_muonMultiplicity(tr.getVec<int>("Jets_muonMultiplicity"))
            , Jets_ID(tr.getVec<bool>("Jets_ID"))
            , JetID(tr.getVar<bool>("JetID"))
            , Jets_partonFlavor(tr.getVec<int>("Jets_partonFlavor"))
            , Jets_ptD(tr.getVec<float>("Jets_ptD"))
            , Jets_axismajor(tr.getVec<float>("Jets_axismajor"))    
            , Jets_axisminor(tr.getVec<float>("Jets_axisminor"))
            , Jets_multiplicity(tr.getVec<int>("Jets_multiplicity"))
            , Jets_neutralEmEnergyFraction(tr.getVec<float>("Jets_neutralEmEnergyFraction"))
            , Jets_chargedEmEnergyFraction(tr.getVec<float>("Jets_chargedEmEnergyFraction"))
            , Jets_neutralHadronEnergyFraction(tr.getVec<float>("Jets_neutralHadronEnergyFraction"))
            , Jets_chargedHadronEnergyFraction(tr.getVec<float>("Jets_chargedHadronEnergyFraction"))
            , Jets_pileupId(tr.getVec<float>("Jets_pileupId"))
        {
        }
    };

    class JetAK8Collection
    {
    public:
        const std::vector<utility::LorentzVector>&              JetsAK8;
        const std::vector<TLorentzVector>&                      JetsAK8_TLV;
        const std::vector<float>&                               JetsAK8_NsubjettinessTau1;
        const std::vector<float>&                               JetsAK8_NsubjettinessTau2;
        const std::vector<float>&                               JetsAK8_NsubjettinessTau3;
        const std::vector<float>&                               JetsAK8_softDropMass;
        const std::vector<float>&                               JetsAK8_axismajor;
        const std::vector<float>&                               JetsAK8_axisminor;
        const std::vector<std::vector<utility::LorentzVector>>& JetsAK8_subjetsNested_LV;
        const std::vector<std::vector<TLorentzVector>>&         JetsAK8_subjetsNested_TLV;
        const std::vector<float>&                               JetsAK8_DeepTagTvsQCD;
        const std::vector<float>&                               JetsAK8_DeepTagWvsQCD;
        const std::vector<float>&                               JetsAK8_DeepTagHbbvsQCD;
        const std::vector<int>&                                 JetsAK8_multiplicity;
    
    JetAK8Collection(const NTupleReader& tr) 
            : JetsAK8(tr.getVec<utility::LorentzVector>("JetsAK8"))
            , JetsAK8_TLV(tr.getVec<TLorentzVector>("JetsAK8_TLV"))
            , JetsAK8_NsubjettinessTau1(tr.getVec<float>("JetsAK8_NsubjettinessTau1"))
            , JetsAK8_NsubjettinessTau2(tr.getVec<float>("JetsAK8_NsubjettinessTau2"))
            , JetsAK8_NsubjettinessTau3(tr.getVec<float>("JetsAK8_NsubjettinessTau3"))
            , JetsAK8_softDropMass(tr.getVec<float>("JetsAK8_softDropMass"))
            , JetsAK8_axismajor(tr.getVec<float>("JetsAK8_axismajor"))
            , JetsAK8_axisminor(tr.getVec<float>("JetsAK8_axisminor"))
            , JetsAK8_subjetsNested_LV(tr.getVec<std::vector<utility::LorentzVector>>("JetsAK8_subjetsNested_LV"))
            , JetsAK8_subjetsNested_TLV(tr.getVec<std::vector<TLorentzVector>>("JetsAK8_subjetsNested_TLV"))
            , JetsAK8_DeepTagTvsQCD(tr.getVec<float>("JetsAK8_DeepTagTvsQCD"))
            , JetsAK8_DeepTagWvsQCD(tr.getVec<float>("JetsAK8_DeepTagWvsQCD"))
            , JetsAK8_DeepTagHbbvsQCD(tr.getVec<float>("JetsAK8_DeepTagHbbvsQCD"))
            , JetsAK8_multiplicity(tr.getVec<int>("JetsAK8_multiplicity"))
        {
        }
    };

    class Factor
    {
    private:
        const std::vector<float>& Jets_jerFactor;
        const std::vector<float>& Jets_jecUnc;
        const std::vector<float>& JetsJECup_jerFactor;
        const std::vector<float>& JetsJECdown_jerFactor;
        const std::vector<float>& Jets_jerFactorUp;
        const std::vector<float>& Jets_jerFactorDown;

    public:
        float factor(const std::string& name, const int i, const int j) const
        {            
            if     (name == "JECup")   return (1./Jets_jerFactor.at(i))*(1+Jets_jecUnc.at(i))*JetsJECup_jerFactor.at(j);
            else if(name == "JECdown") return (1./Jets_jerFactor.at(i))*(1-Jets_jecUnc.at(i))*JetsJECdown_jerFactor.at(j);
            else if(name == "JERup")   return (1./Jets_jerFactor.at(i))*Jets_jerFactorUp.at(i);
            else if(name == "JERdown") return (1./Jets_jerFactor.at(i))*Jets_jerFactorDown.at(i);                
            else                       return -9999.9;
        }
        
        Factor(const NTupleReader& tr, const std::string& name = "") 
            : Jets_jerFactor(tr.getVec<float>("Jets"+name+"_jerFactor"))
            , Jets_jecUnc(tr.getVec<float>("Jets"+name+"_jecUnc"))
            , JetsJECup_jerFactor(tr.getVec<float>("Jets"+name+"JECup_jerFactor"))
            , JetsJECdown_jerFactor(tr.getVec<float>("Jets"+name+"JECdown_jerFactor"))
            , Jets_jerFactorUp(tr.getVec<float>("Jets"+name+"_jerFactorUp"))
            , Jets_jerFactorDown(tr.getVec<float>("Jets"+name+"_jerFactorDown"))
        {
        }
    };

    void deriveJetCollection(NTupleReader& tr, const JetCollection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("Jets"+name+"_origIndex");
               
        auto& newJets                             = tr.createDerivedVec<utility::LorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_TLV                         = tr.createDerivedVec<TLorentzVector>("Jets_TLV"+name, jc.Jets.size());
        auto& newJets_bDiscriminatorCSV           = tr.createDerivedVec<float>("Jets"+name+"_bDiscriminatorCSV", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb         = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb        = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb          = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb    = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourproblepb  = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourproblepb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds   = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobc         = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobudsg      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobudsg", jc.Jets.size());
        auto& newJets_qgLikelihood                = tr.createDerivedVec<float>("Jets"+name+"_qgLikelihood", jc.Jets.size());
        auto& newJets_ID                          = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor                = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());
        auto& newJets_ptD                         = tr.createDerivedVec<float>("Jets"+name+"_ptD", jc.Jets.size());
        auto& newJets_axismajor                   = tr.createDerivedVec<float>("Jets"+name+"_axismajor", jc.Jets.size());
        auto& newJets_axisminor                   = tr.createDerivedVec<float>("Jets"+name+"_axisminor", jc.Jets.size());
        auto& newJets_multiplicity                = tr.createDerivedVec<int>("Jets"+name+"_multiplicity", jc.Jets.size());
        auto& newJets_neutralEmEnergyFraction     = tr.createDerivedVec<float>("Jets"+name+"_neutralEmEnergyFraction", jc.Jets.size());
        auto& newJets_chargedEmEnergyFraction     = tr.createDerivedVec<float>("Jets"+name+"_chargedEmEnergyFraction", jc.Jets.size());
        auto& newJets_neutralHadronEnergyFraction = tr.createDerivedVec<float>("Jets"+name+"_neutralHadronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronEnergyFraction = tr.createDerivedVec<float>("Jets"+name+"_chargedHadronEnergyFraction", jc.Jets.size());
        auto& newJets_muonEnergyFraction          = tr.createDerivedVec<float>("Jets"+name+"_muonEnergyFraction", jc.Jets.size());
        auto& newJets_hfHadronEnergyFraction      = tr.createDerivedVec<float>("Jets"+name+"_hfHadronEnergyFraction", jc.Jets.size());
        auto& newJets_hfEMEnergyFraction          = tr.createDerivedVec<float>("Jets"+name+"_hfEMEnergyFraction", jc.Jets.size());
        auto& newJets_photonEnergyFraction        = tr.createDerivedVec<float>("Jets"+name+"_photonEnergyFraction", jc.Jets.size());
        auto& newJets_electronEnergyFraction      = tr.createDerivedVec<float>("Jets"+name+"_electronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_chargedHadronMultiplicity", jc.Jets.size());
        auto& newJets_neutralHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_neutralHadronMultiplicity", jc.Jets.size());
        auto& newJets_photonMultiplicity          = tr.createDerivedVec<int>("Jets"+name+"_photonMultiplicity", jc.Jets.size());
        auto& newJets_electronMultiplicity        = tr.createDerivedVec<int>("Jets"+name+"_electronMultiplicity", jc.Jets.size());
        auto& newJets_muonMultiplicity            = tr.createDerivedVec<int>("Jets"+name+"_muonMultiplicity", jc.Jets.size());
        auto& newJets_pileupId                    = tr.createDerivedVec<float>("Jets"+name+"_pileupId", jc.Jets.size());

        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i    = newIndex[newJets_origIndex.at(j)];
            auto lv  = jc.Jets.at(i)*f.factor(name,i,j);
            auto tlv = jc.Jets_TLV.at(i)*f.factor(name,i,j);

            newJets.at(j).SetPt(lv.Pt()); newJets.at(j).SetEta(lv.Eta()); newJets.at(j).SetPhi(lv.Phi()); newJets.at(j).SetE(lv.E());
            newJets_TLV.at(j).SetPtEtaPhiE(tlv.Pt(), tlv.Eta(), tlv.Phi(), tlv.E());
            newJets_bDiscriminatorCSV.at(j)           = jc.Jets_bDiscriminatorCSV.at(i);
            newJets_bJetTagDeepCSVprobb.at(j)         = jc.Jets_bJetTagDeepCSVprobb.at(i);
            newJets_bJetTagDeepCSVprobbb.at(j)        = jc.Jets_bJetTagDeepCSVprobbb.at(i);
            newJets_bJetTagDeepCSVtotb.at(j)          = ( jc.Jets_bJetTagDeepCSVprobb.at(i) + jc.Jets_bJetTagDeepCSVprobbb.at(i) );
            newJets_bJetTagDeepFlavourprobb.at(j)     = jc.Jets_bJetTagDeepFlavourprobb.at(i);
            newJets_bJetTagDeepFlavourprobbb.at(j)    = jc.Jets_bJetTagDeepFlavourprobbb.at(i);
            newJets_bJetTagDeepFlavourproblepb.at(j)  = jc.Jets_bJetTagDeepFlavourproblepb.at(i);
            newJets_bJetTagDeepFlavourtotb.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobb.at(i) + jc.Jets_bJetTagDeepFlavourprobbb.at(i) + jc.Jets_bJetTagDeepFlavourproblepb.at(i) );
            newJets_bJetTagDeepFlavourtotq.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobc.at(i) + jc.Jets_bJetTagDeepFlavourprobuds.at(i) );
            newJets_bJetTagDeepFlavourprobg.at(j)     = jc.Jets_bJetTagDeepFlavourprobg.at(i); 
            newJets_bJetTagDeepFlavourprobc.at(j)     = jc.Jets_bJetTagDeepFlavourprobc.at(i); 
            newJets_bJetTagDeepFlavourprobuds.at(j)   = jc.Jets_bJetTagDeepFlavourprobuds.at(i); 
            newJets_bJetTagDeepCSVprobc.at(j)         = jc.Jets_bJetTagDeepCSVprobc.at(i);
            newJets_bJetTagDeepCSVprobudsg.at(j)      = jc.Jets_bJetTagDeepCSVprobudsg.at(i);
            newJets_muonEnergyFraction.at(j)          = jc.Jets_muonEnergyFraction.at(i);
            newJets_hfHadronEnergyFraction.at(j)      = jc.Jets_hfHadronEnergyFraction.at(i);
            newJets_hfEMEnergyFraction.at(j)          = jc.Jets_hfEMEnergyFraction.at(i);
            newJets_photonEnergyFraction.at(j)        = jc.Jets_photonEnergyFraction.at(i);
            newJets_electronEnergyFraction.at(j)      = jc.Jets_electronEnergyFraction.at(i);
            newJets_chargedHadronMultiplicity.at(j)   = jc.Jets_chargedHadronMultiplicity.at(i);
            newJets_neutralHadronMultiplicity.at(j)   = jc.Jets_neutralHadronMultiplicity.at(i);
            newJets_photonMultiplicity.at(j)          = jc.Jets_photonMultiplicity.at(i);
            newJets_electronMultiplicity.at(j)        = jc.Jets_electronMultiplicity.at(i);
            newJets_muonMultiplicity.at(j)            = jc.Jets_muonMultiplicity.at(i);
            newJets_qgLikelihood.at(j)                = jc.Jets_qgLikelihood.at(i); 
            newJets_ID.at(j)                          = jc.Jets_ID.at(i);
            newJets_partonFlavor.at(j)                = jc.Jets_partonFlavor.at(i);
            newJets_ptD.at(j)                         = jc.Jets_ptD.at(i);
            newJets_axismajor.at(j)                   = jc.Jets_axismajor.at(i);
            newJets_axisminor.at(j)                   = jc.Jets_axisminor.at(i);
            newJets_multiplicity.at(j)                = jc.Jets_multiplicity.at(i);
            newJets_neutralEmEnergyFraction.at(j)     = jc.Jets_neutralEmEnergyFraction.at(i);
            newJets_chargedEmEnergyFraction.at(j)     = jc.Jets_chargedEmEnergyFraction.at(i);
            newJets_neutralHadronEnergyFraction.at(j) = jc.Jets_neutralHadronEnergyFraction.at(i);
            newJets_chargedHadronEnergyFraction.at(j) = jc.Jets_chargedHadronEnergyFraction.at(i);
            if (isinf(jc.Jets_pileupId.at(i))){
                newJets_pileupId.at(j)                    = 999.0;
            } else {
                newJets_pileupId.at(j)                    = jc.Jets_pileupId.at(i);
            }
        }
    }

    void derivePtMassScaledJetCollection(NTupleReader& tr, const JetCollection& jc, const std::string& name, double scalePt, double scaleMass)
    {
        tr.registerDerivedVar("JetID"+name, jc.JetID);        

        auto& newJets                             = tr.createDerivedVec<utility::LorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_TLV                         = tr.createDerivedVec<TLorentzVector>("Jets_TLV"+name, jc.Jets.size());
        auto& newJets_bDiscriminatorCSV           = tr.createDerivedVec<float>("Jets"+name+"_bDiscriminatorCSV", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb         = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb        = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb          = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb    = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourproblepb  = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourproblepb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc     = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds   = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());
        auto& newJets_qgLikelihood                = tr.createDerivedVec<float>("Jets"+name+"_qgLikelihood", jc.Jets.size());
        auto& newJets_ID                          = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor                = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());
        auto& newJets_ptD                         = tr.createDerivedVec<float>("Jets"+name+"_ptD", jc.Jets.size());
        auto& newJets_axismajor                   = tr.createDerivedVec<float>("Jets"+name+"_axismajor", jc.Jets.size());
        auto& newJets_axisminor                   = tr.createDerivedVec<float>("Jets"+name+"_axisminor", jc.Jets.size());
        auto& newJets_multiplicity                = tr.createDerivedVec<int>("Jets"+name+"_multiplicity", jc.Jets.size());
        auto& newJets_neutralEmEnergyFraction     = tr.createDerivedVec<float>("Jets"+name+"_neutralEmEnergyFraction", jc.Jets.size());
        auto& newJets_chargedEmEnergyFraction     = tr.createDerivedVec<float>("Jets"+name+"_chargedEmEnergyFraction", jc.Jets.size());
        auto& newJets_neutralHadronEnergyFraction = tr.createDerivedVec<float>("Jets"+name+"_neutralHadronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronEnergyFraction = tr.createDerivedVec<float>("Jets"+name+"_chargedHadronEnergyFraction", jc.Jets.size());
        auto& newJets_muonEnergyFraction          = tr.createDerivedVec<float>("Jets"+name+"_muonEnergyFraction", jc.Jets.size());
        auto& newJets_hfHadronEnergyFraction      = tr.createDerivedVec<float>("Jets"+name+"_hfHadronEnergyFraction", jc.Jets.size());
        auto& newJets_hfEMEnergyFraction          = tr.createDerivedVec<float>("Jets"+name+"_hfEMEnergyFraction", jc.Jets.size());
        auto& newJets_photonEnergyFraction        = tr.createDerivedVec<float>("Jets"+name+"_photonEnergyFraction", jc.Jets.size());
        auto& newJets_electronEnergyFraction      = tr.createDerivedVec<float>("Jets"+name+"_electronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_chargedHadronMultiplicity", jc.Jets.size());
        auto& newJets_neutralHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_neutralHadronMultiplicity", jc.Jets.size());
        auto& newJets_photonMultiplicity          = tr.createDerivedVec<int>("Jets"+name+"_photonMultiplicity", jc.Jets.size());
        auto& newJets_electronMultiplicity        = tr.createDerivedVec<int>("Jets"+name+"_electronMultiplicity", jc.Jets.size());
        auto& newJets_muonMultiplicity            = tr.createDerivedVec<int>("Jets"+name+"_muonMultiplicity", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobc         = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobudsg      = tr.createDerivedVec<float>("Jets"+name+"_bJetTagDeepCSVprobudsg", jc.Jets.size());
        auto& newJets_pileupId                    = tr.createDerivedVec<float>("Jets"+name+"_pileupId", jc.Jets.size());

        for(unsigned j = 0; j < jc.Jets.size(); ++j)
        {
            newJets.at(j).SetPt( scalePt*jc.Jets[j].Pt()); newJets.at(j).SetEta(jc.Jets[j].Eta()); newJets.at(j).SetPhi(jc.Jets[j].Phi()); newJets.at(j).SetE(scaleMass*jc.Jets[j].E());
            newJets_TLV.at(j).SetPtEtaPhiE( scalePt*jc.Jets_TLV[j].Pt(), jc.Jets_TLV[j].Eta(), jc.Jets_TLV[j].Phi(), scaleMass*jc.Jets_TLV[j].E());
            newJets_bDiscriminatorCSV.at(j)           = jc.Jets_bDiscriminatorCSV.at(j);
            newJets_bJetTagDeepCSVprobb.at(j)         = jc.Jets_bJetTagDeepCSVprobb.at(j);
            newJets_bJetTagDeepCSVprobbb.at(j)        = jc.Jets_bJetTagDeepCSVprobbb.at(j);
            newJets_bJetTagDeepCSVtotb.at(j)          = ( jc.Jets_bJetTagDeepCSVprobb.at(j) + jc.Jets_bJetTagDeepCSVprobbb.at(j) );
            newJets_bJetTagDeepFlavourprobb.at(j)     = jc.Jets_bJetTagDeepFlavourprobb.at(j);
            newJets_bJetTagDeepFlavourprobbb.at(j)    = jc.Jets_bJetTagDeepFlavourprobbb.at(j);
            newJets_bJetTagDeepFlavourproblepb.at(j)  = jc.Jets_bJetTagDeepFlavourproblepb.at(j);
            newJets_bJetTagDeepFlavourtotb.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobb.at(j) + jc.Jets_bJetTagDeepFlavourprobbb.at(j) + jc.Jets_bJetTagDeepFlavourproblepb.at(j) );
            newJets_bJetTagDeepFlavourtotq.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobc.at(j) + jc.Jets_bJetTagDeepFlavourprobuds.at(j) );
            newJets_bJetTagDeepFlavourprobg.at(j)     = jc.Jets_bJetTagDeepFlavourprobg.at(j); 
            newJets_bJetTagDeepFlavourprobc.at(j)     = jc.Jets_bJetTagDeepFlavourprobc.at(j); 
            newJets_bJetTagDeepFlavourprobuds.at(j)   = jc.Jets_bJetTagDeepFlavourprobuds.at(j); 
            newJets_muonEnergyFraction.at(j)          = jc.Jets_muonEnergyFraction.at(j);
            newJets_hfHadronEnergyFraction.at(j)      = jc.Jets_hfHadronEnergyFraction.at(j);
            newJets_hfEMEnergyFraction.at(j)          = jc.Jets_hfEMEnergyFraction.at(j);
            newJets_photonEnergyFraction.at(j)        = jc.Jets_photonEnergyFraction.at(j);
            newJets_electronEnergyFraction.at(j)      = jc.Jets_electronEnergyFraction.at(j);
            newJets_chargedHadronMultiplicity.at(j)   = jc.Jets_chargedHadronMultiplicity.at(j);
            newJets_neutralHadronMultiplicity.at(j)   = jc.Jets_neutralHadronMultiplicity.at(j);
            newJets_photonMultiplicity.at(j)          = jc.Jets_photonMultiplicity.at(j);
            newJets_electronMultiplicity.at(j)        = jc.Jets_electronMultiplicity.at(j);
            newJets_muonMultiplicity.at(j)            = jc.Jets_muonMultiplicity.at(j);
            newJets_bJetTagDeepCSVprobc.at(j)         = jc.Jets_bJetTagDeepCSVprobc.at(j);
            newJets_bJetTagDeepCSVprobudsg.at(j)      = jc.Jets_bJetTagDeepCSVprobudsg.at(j);
            newJets_qgLikelihood.at(j)                = jc.Jets_qgLikelihood.at(j); 
            newJets_ID.at(j)                          = jc.Jets_ID.at(j);
            newJets_partonFlavor.at(j)                = jc.Jets_partonFlavor.at(j);
            newJets_ptD.at(j)                         = jc.Jets_ptD.at(j);
            newJets_axismajor.at(j)                   = jc.Jets_axismajor.at(j);
            newJets_axisminor.at(j)                   = jc.Jets_axisminor.at(j);
            newJets_multiplicity.at(j)                = jc.Jets_multiplicity.at(j);
            newJets_neutralEmEnergyFraction.at(j)     = jc.Jets_neutralEmEnergyFraction.at(j);
            newJets_chargedEmEnergyFraction.at(j)     = jc.Jets_chargedEmEnergyFraction.at(j);
            newJets_neutralHadronEnergyFraction.at(j) = jc.Jets_neutralHadronEnergyFraction.at(j);
            newJets_chargedHadronEnergyFraction.at(j) = jc.Jets_chargedHadronEnergyFraction.at(j);
            if (isinf(jc.Jets_pileupId.at(j))){
                newJets_pileupId.at(j)                    = 999.0;
            } else {
                newJets_pileupId.at(j)                    = jc.Jets_pileupId.at(j);
            }
        }
    }

    void deriveJetAK8Collection(NTupleReader& tr, const JetAK8Collection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("JetsAK8"+name+"_origIndex");

        auto& newJetsAK8                    = tr.createDerivedVec<utility::LorentzVector>("JetsAK8"+name, jc.JetsAK8.size());
        auto& newJetsAK8_TLV                = tr.createDerivedVec<TLorentzVector>("JetsAK8_TLV"+name, jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau1  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau1", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau2  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau2", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau3  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau3", jc.JetsAK8.size());
        auto& newJetsAK8_softDropMass       = tr.createDerivedVec<float>("JetsAK8"+name+"_softDropMass", jc.JetsAK8.size());
        auto& newJetsAK8_axismajor          = tr.createDerivedVec<float>("JetsAK8"+name+"_axismajor", jc.JetsAK8.size());
        auto& newJetsAK8_axisminor          = tr.createDerivedVec<float>("JetsAK8"+name+"_axisminor", jc.JetsAK8.size());
        auto& newJetsAK8_subjetsNested_LV   = tr.createDerivedVec<std::vector<utility::LorentzVector>>("JetsAK8"+name+"_subjetsNested_LV", jc.JetsAK8.size());
        auto& newJetsAK8_subjetsNested_TLV  = tr.createDerivedVec<std::vector<TLorentzVector>>("JetsAK8"+name+"_subjetsNested_TLV", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagTvsQCD      = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagTvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagWvsQCD      = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagWvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagHbbvsQCD    = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagHbbvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_multiplicity       = tr.createDerivedVec<int>("JetsAK8"+name+"_multiplicity", jc.JetsAK8.size());
             
        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i    = newIndex[newJets_origIndex.at(j)];
            auto lv  = jc.JetsAK8.at(i)*f.factor(name,i,j);
            auto tlv = jc.JetsAK8_TLV.at(i)*f.factor(name,i,j);

            newJetsAK8.at(j).SetPt(lv.Pt()); newJetsAK8.at(j).SetEta(lv.Eta()); newJetsAK8.at(j).SetPhi(lv.Phi()); newJetsAK8.at(j).SetE(lv.E());
            newJetsAK8_TLV.at(j).SetPtEtaPhiE(tlv.Pt(), tlv.Eta(), tlv.Phi(), tlv.E());
            newJetsAK8_NsubjettinessTau1.at(j)  = jc.JetsAK8_NsubjettinessTau1.at(i);
            newJetsAK8_NsubjettinessTau2.at(j)  = jc.JetsAK8_NsubjettinessTau2.at(i);
            newJetsAK8_NsubjettinessTau3.at(j)  = jc.JetsAK8_NsubjettinessTau3.at(i);
            newJetsAK8_softDropMass.at(j)       = jc.JetsAK8_softDropMass.at(i);
            newJetsAK8_axismajor.at(j)          = jc.JetsAK8_axismajor.at(i);
            newJetsAK8_axisminor.at(j)          = jc.JetsAK8_axisminor.at(i);
            newJetsAK8_subjetsNested_LV.at(j)   = jc.JetsAK8_subjetsNested_LV.at(i);
            newJetsAK8_subjetsNested_TLV.at(j)  = jc.JetsAK8_subjetsNested_TLV.at(i);
            newJetsAK8_DeepTagTvsQCD.at(j)      = jc.JetsAK8_DeepTagTvsQCD.at(i);
            newJetsAK8_DeepTagWvsQCD.at(j)      = jc.JetsAK8_DeepTagWvsQCD.at(i);
            newJetsAK8_DeepTagHbbvsQCD.at(j)    = jc.JetsAK8_DeepTagHbbvsQCD.at(i);
            newJetsAK8_multiplicity.at(j)       = jc.JetsAK8_multiplicity.at(i);
        }
    }

    void derivePtMassScaledJetAK8Collection(NTupleReader& tr, const JetAK8Collection& jc, const std::string& name, double scalePt = 1.0, double scaleMass = 1.0)
    {

        auto& newJetsAK8                    = tr.createDerivedVec<utility::LorentzVector>("JetsAK8"+name, jc.JetsAK8.size());
        auto& newJetsAK8_TLV                = tr.createDerivedVec<TLorentzVector>("JetsAK8_TLV"+name, jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau1  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau1", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau2  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau2", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau3  = tr.createDerivedVec<float>("JetsAK8"+name+"_NsubjettinessTau3", jc.JetsAK8.size());
        auto& newJetsAK8_softDropMass       = tr.createDerivedVec<float>("JetsAK8"+name+"_softDropMass", jc.JetsAK8.size());
        auto& newJetsAK8_axismajor          = tr.createDerivedVec<float>("JetsAK8"+name+"_axismajor", jc.JetsAK8.size());
        auto& newJetsAK8_axisminor          = tr.createDerivedVec<float>("JetsAK8"+name+"_axisminor", jc.JetsAK8.size());
        auto& newJetsAK8_subjetsNested_LV   = tr.createDerivedVec<std::vector<utility::LorentzVector>>("JetsAK8"+name+"_subjetsNested_LV", jc.JetsAK8.size());
        auto& newJetsAK8_subjetsNested_TLV  = tr.createDerivedVec<std::vector<TLorentzVector>>("JetsAK8"+name+"_subjetsNested_TLV", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagTvsQCD      = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagTvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagWvsQCD      = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagWvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_DeepTagHbbvsQCD    = tr.createDerivedVec<float>("JetsAK8"+name+"_DeepTagHbbvsQCD", jc.JetsAK8.size());
        auto& newJetsAK8_multiplicity       = tr.createDerivedVec<int>("JetsAK8"+name+"_multiplicity", jc.JetsAK8.size());
             
        for(unsigned j = 0; j < jc.JetsAK8.size(); ++j)
        {
            newJetsAK8.at(j).SetPt( scalePt*jc.JetsAK8[j].Pt()); newJetsAK8.at(j).SetEta(jc.JetsAK8[j].Eta()); newJetsAK8.at(j).SetPhi(jc.JetsAK8[j].Phi()); newJetsAK8.at(j).SetE( scaleMass*jc.JetsAK8[j].E());
            newJetsAK8_TLV.at(j).SetPtEtaPhiE( scalePt*jc.JetsAK8_TLV[j].Pt(), jc.JetsAK8_TLV[j].Eta(), jc.JetsAK8_TLV[j].Phi(), scaleMass*jc.JetsAK8_TLV[j].E());

            newJetsAK8_NsubjettinessTau1.at(j)  = jc.JetsAK8_NsubjettinessTau1.at(j);
            newJetsAK8_NsubjettinessTau2.at(j)  = jc.JetsAK8_NsubjettinessTau2.at(j);
            newJetsAK8_NsubjettinessTau3.at(j)  = jc.JetsAK8_NsubjettinessTau3.at(j);
            newJetsAK8_softDropMass.at(j)       = jc.JetsAK8_softDropMass.at(j);
            newJetsAK8_axismajor.at(j)          = jc.JetsAK8_axismajor.at(j);
            newJetsAK8_axisminor.at(j)          = jc.JetsAK8_axisminor.at(j);
            newJetsAK8_subjetsNested_LV.at(j)   = jc.JetsAK8_subjetsNested_LV.at(j);
            newJetsAK8_subjetsNested_TLV.at(j)  = jc.JetsAK8_subjetsNested_TLV.at(j);
            newJetsAK8_DeepTagTvsQCD.at(j)      = jc.JetsAK8_DeepTagTvsQCD.at(j);
            newJetsAK8_DeepTagWvsQCD.at(j)      = jc.JetsAK8_DeepTagWvsQCD.at(j);
            newJetsAK8_DeepTagHbbvsQCD.at(j)    = jc.JetsAK8_DeepTagHbbvsQCD.at(j);
            newJetsAK8_multiplicity.at(j)       = jc.JetsAK8_multiplicity.at(j);
        }
    }

    void prepNTupleVars(NTupleReader& tr)
    {
        // Make nested vector of vectors for AK8 subjets
        const auto& JetsAK8_subjetsCounts = tr.getVec<int>("JetsAK8_subjetsCounts");
        const auto& JetsAK8_subjets       = tr.getVec<utility::LorentzVector>("JetsAK8_subjets");

        const auto& Jets    = tr.getVec<utility::LorentzVector>("Jets");
        const auto& JetsAK8 = tr.getVec<utility::LorentzVector>("JetsAK8");

        auto& JetsAK8_subjetsNested_LV  = tr.createDerivedVec<std::vector<utility::LorentzVector>>("JetsAK8_subjetsNested_LV", JetsAK8_subjetsCounts.size());
        auto& JetsAK8_subjetsNested_TLV = tr.createDerivedVec<std::vector<TLorentzVector>>("JetsAK8_subjetsNested_TLV", JetsAK8_subjetsCounts.size());

        JetsAK8_subjetsNested_LV  = utility::nestVecOfVec<utility::LorentzVector, utility::LorentzVector>(JetsAK8_subjets, JetsAK8_subjetsCounts);
        JetsAK8_subjetsNested_TLV = utility::nestVecOfVec<TLorentzVector, utility::LorentzVector>(JetsAK8_subjets, JetsAK8_subjetsCounts);

        auto& JetsTLV = tr.createDerivedVec<TLorentzVector>("Jets_TLV", Jets.size());
        JetsTLV = utility::convertVectorOfLV<TLorentzVector, utility::LorentzVector>(Jets);
    
        auto& JetsAK8TLV = tr.createDerivedVec<TLorentzVector>("JetsAK8_TLV", JetsAK8.size());
        JetsAK8TLV = utility::convertVectorOfLV<TLorentzVector, utility::LorentzVector>(JetsAK8);

        // Creating the jet pT and mass scaled collection
        JetCollection jc(tr);
        JetAK8Collection jcAK8(tr);
        const auto& runYear = tr.getVar<std::string>("runYear");
        derivePtMassScaledJetCollection(tr, jc, "mpTScaled", pTMass_[runYear].first, pTMass_[runYear].second);
        derivePtMassScaledJetAK8Collection(tr, jcAK8, "mpTScaled", 1.0, 1.0);

        // Create DeepCSV b-jet discriminator vector
        const auto& Jets_bJetTagDeepCSVprobb        = tr.getVec<float>("Jets_bJetTagDeepCSVprobb");
        const auto& Jets_bJetTagDeepCSVprobbb       = tr.getVec<float>("Jets_bJetTagDeepCSVprobbb");
              auto& Jets_bJetTagDeepCSVtotb         = tr.createDerivedVec<float>("Jets_bJetTagDeepCSVtotb", Jets_bJetTagDeepCSVprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobb    = tr.getVec<float>("Jets_bJetTagDeepFlavourprobb");
        const auto& Jets_bJetTagDeepFlavourprobbb   = tr.getVec<float>("Jets_bJetTagDeepFlavourprobbb");
        const auto& Jets_bJetTagDeepFlavourproblepb = tr.getVec<float>("Jets_bJetTagDeepFlavourproblepb");
              auto& Jets_bJetTagDeepFlavourtotb     = tr.createDerivedVec<float>("Jets_bJetTagDeepFlavourtotb", Jets_bJetTagDeepFlavourprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobc    = tr.getVec<float>("Jets_bJetTagDeepFlavourprobc");
        const auto& Jets_bJetTagDeepFlavourprobuds  = tr.getVec<float>("Jets_bJetTagDeepFlavourprobuds");
              auto& Jets_bJetTagDeepFlavourtotq     = tr.createDerivedVec<float>("Jets_bJetTagDeepFlavourtotq", Jets_bJetTagDeepFlavourprobc.size());

        for(unsigned j = 0; j < Jets_bJetTagDeepCSVprobb.size(); ++j)
        {
            Jets_bJetTagDeepCSVtotb.at(j)     = ( Jets_bJetTagDeepCSVprobb.at(j) + Jets_bJetTagDeepCSVprobbb.at(j) );
            Jets_bJetTagDeepFlavourtotb.at(j) = ( Jets_bJetTagDeepFlavourprobb.at(j) + Jets_bJetTagDeepFlavourprobbb.at(j) + Jets_bJetTagDeepFlavourproblepb.at(j) );
            Jets_bJetTagDeepFlavourtotq.at(j) = ( Jets_bJetTagDeepFlavourprobc.at(j) + Jets_bJetTagDeepFlavourprobuds.at(j) );
        }
 
        // Create JEC/R variation variables if needed
        const auto& runtype = tr.getVar<std::string>("runtype");
        if( !tr.checkBranchInTree("JetsJECup") && runtype == "MC")
        {
            //Make AK4 jet JEC/R variables
            const auto& Jets_origIndex = tr.getVec<int>("Jets_origIndex");
            std::vector<int> newIndex(Jets_origIndex.size());
            for(unsigned j = 0; j < Jets_origIndex.size(); ++j)
            {
                newIndex[Jets_origIndex.at(j)] = j;
            }
            Factor f(tr);
            deriveJetCollection(tr, jc, f, newIndex, "JECup");
            deriveJetCollection(tr, jc, f, newIndex, "JECdown");
            deriveJetCollection(tr, jc, f, newIndex, "JERup");
            deriveJetCollection(tr, jc, f, newIndex, "JERdown");

            //Make AK8 jet JEC/R variables
            const auto& JetsAK8_origIndex = tr.getVec<int>("JetsAK8_origIndex");
            std::vector<int> newIndexAK8(JetsAK8_origIndex.size());
            for(unsigned j = 0; j < JetsAK8_origIndex.size(); ++j)
            {
                newIndexAK8[JetsAK8_origIndex.at(j)] = j;
            }
            Factor fAK8(tr, "AK8");
            deriveJetAK8Collection(tr, jcAK8, fAK8, newIndexAK8, "JECup");
            deriveJetAK8Collection(tr, jcAK8, fAK8, newIndexAK8, "JECdown");
            deriveJetAK8Collection(tr, jcAK8, fAK8, newIndexAK8, "JERup");
            deriveJetAK8Collection(tr, jcAK8, fAK8, newIndexAK8, "JERdown");
        }        
    }
public:
    //PrepNTupleVars() : pTMass_({{"2016preVFP",{0.95,0.95}},{"2016postVFP",{0.95,0.95}},{"2017",{0.95,1.01}},{"2018pre",{0.95,0.98}},{"2018post",{0.95,0.98}}})
    PrepNTupleVars() : pTMass_({{"2016preVFP",{1,1}},{"2016postVFP",{1,1}},{"2017",{1,1}},{"2018pre",{1,1}},{"2018post",{1,1}}, {"2018", {1,1}}})
    {
    }

    void operator()(NTupleReader& tr)
    {
        prepNTupleVars(tr);
    }
};

#endif
