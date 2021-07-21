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
        const std::vector<TLorentzVector>& Jets;
        const std::vector<double>&         Jets_bDiscriminatorCSV;
        const std::vector<double>&         Jets_bJetTagDeepCSVprobb;
        const std::vector<double>&         Jets_bJetTagDeepCSVprobbb;
        const std::vector<double>&         Jets_bJetTagDeepFlavourprobb;
        const std::vector<double>&         Jets_bJetTagDeepFlavourprobbb;
        const std::vector<double>&         Jets_bJetTagDeepFlavourprobc;
        const std::vector<double>&         Jets_bJetTagDeepFlavourprobuds;
        const std::vector<double>&         Jets_bJetTagDeepFlavourprobg;
        const std::vector<double>&         Jets_bJetTagDeepCSVprobc;
        const std::vector<double>&         Jets_bJetTagDeepCSVprobudsg;
        const std::vector<double>&         Jets_qgLikelihood;
        const std::vector<double>&         Jets_muonEnergyFraction;
        const std::vector<double>&         Jets_hfHadronEnergyFraction;
        const std::vector<double>&         Jets_hfEMEnergyFraction;
        const std::vector<double>&         Jets_photonEnergyFraction;
        const std::vector<double>&         Jets_electronEnergyFraction;
        const std::vector<int>&            Jets_chargedHadronMultiplicity;
        const std::vector<int>&            Jets_neutralHadronMultiplicity;
        const std::vector<int>&            Jets_photonMultiplicity;
        const std::vector<int>&            Jets_electronMultiplicity;
        const std::vector<int>&            Jets_muonMultiplicity;
        const std::vector<bool>&           Jets_ID;
        const bool&                        JetID;
        const std::vector<int>&            Jets_partonFlavor;
        const std::vector<double>&         Jets_ptD;
        const std::vector<double>&         Jets_axismajor;
        const std::vector<double>&         Jets_axisminor;
        const std::vector<int>&            Jets_multiplicity;
        const std::vector<double>&         Jets_neutralEmEnergyFraction;
        const std::vector<double>&         Jets_chargedEmEnergyFraction;
        const std::vector<double>&         Jets_neutralHadronEnergyFraction;
        const std::vector<double>&         Jets_chargedHadronEnergyFraction;
        
    JetCollection(const NTupleReader& tr) 
            : Jets(tr.getVec<TLorentzVector>("Jets"))
            , Jets_bDiscriminatorCSV(tr.getVec<double>("Jets_bDiscriminatorCSV")) 
            , Jets_bJetTagDeepCSVprobb(tr.getVec<double>("Jets_bJetTagDeepCSVprobb"))
            , Jets_bJetTagDeepCSVprobbb(tr.getVec<double>("Jets_bJetTagDeepCSVprobbb"))
            , Jets_bJetTagDeepFlavourprobb(tr.getVec<double>("Jets_bJetTagDeepFlavourprobb"))
            , Jets_bJetTagDeepFlavourprobbb(tr.getVec<double>("Jets_bJetTagDeepFlavourprobbb"))
            , Jets_bJetTagDeepFlavourprobc(tr.getVec<double>("Jets_bJetTagDeepFlavourprobc"))
            , Jets_bJetTagDeepFlavourprobuds(tr.getVec<double>("Jets_bJetTagDeepFlavourprobuds"))
            , Jets_bJetTagDeepFlavourprobg(tr.getVec<double>("Jets_bJetTagDeepFlavourprobg"))
            , Jets_bJetTagDeepCSVprobc(tr.getVec<double>("Jets_bJetTagDeepCSVprobc"))
            , Jets_bJetTagDeepCSVprobudsg(tr.getVec<double>("Jets_bJetTagDeepCSVprobudsg"))
            , Jets_qgLikelihood(tr.getVec<double>("Jets_qgLikelihood"))
            , Jets_muonEnergyFraction(tr.getVec<double>("Jets_muonEnergyFraction"))
            , Jets_hfHadronEnergyFraction(tr.getVec<double>("Jets_hfHadronEnergyFraction"))
            , Jets_hfEMEnergyFraction(tr.getVec<double>("Jets_hfEMEnergyFraction"))
            , Jets_photonEnergyFraction(tr.getVec<double>("Jets_photonEnergyFraction"))
            , Jets_electronEnergyFraction(tr.getVec<double>("Jets_electronEnergyFraction"))
            , Jets_chargedHadronMultiplicity(tr.getVec<int>("Jets_chargedHadronMultiplicity"))
            , Jets_neutralHadronMultiplicity(tr.getVec<int>("Jets_neutralHadronMultiplicity"))
            , Jets_photonMultiplicity(tr.getVec<int>("Jets_photonMultiplicity"))
            , Jets_electronMultiplicity(tr.getVec<int>("Jets_electronMultiplicity"))
            , Jets_muonMultiplicity(tr.getVec<int>("Jets_muonMultiplicity"))
            , Jets_ID(tr.getVec<bool>("Jets_ID"))
            , JetID(tr.getVar<bool>("JetID"))
            , Jets_partonFlavor(tr.getVec<int>("Jets_partonFlavor"))
            , Jets_ptD(tr.getVec<double>("Jets_ptD"))
            , Jets_axismajor(tr.getVec<double>("Jets_axismajor"))    
            , Jets_axisminor(tr.getVec<double>("Jets_axisminor"))
            , Jets_multiplicity(tr.getVec<int>("Jets_multiplicity"))
            , Jets_neutralEmEnergyFraction(tr.getVec<double>("Jets_neutralEmEnergyFraction"))
            , Jets_chargedEmEnergyFraction(tr.getVec<double>("Jets_chargedEmEnergyFraction"))
            , Jets_neutralHadronEnergyFraction(tr.getVec<double>("Jets_neutralHadronEnergyFraction"))
            , Jets_chargedHadronEnergyFraction(tr.getVec<double>("Jets_chargedHadronEnergyFraction"))
        {
        }
    };

    class JetAK8Collection
    {
    public:
        const std::vector<TLorentzVector>&              JetsAK8;
        const std::vector<double>&                      JetsAK8_NsubjettinessTau1;
        const std::vector<double>&                      JetsAK8_NsubjettinessTau2;
        const std::vector<double>&                      JetsAK8_NsubjettinessTau3;
        const std::vector<double>&                      JetsAK8_softDropMass;
        const std::vector<double>&                      JetsAK8_prunedMass;
        const std::vector<double>&                      JetsAK8_axismajor;
        const std::vector<double>&                      JetsAK8_axisminor;
        const std::vector<std::vector<TLorentzVector>>& JetsAK8_subjets;
        const std::vector<double>&                      JetsAK8_tDiscriminatorDeep;
        const std::vector<double>&                      JetsAK8_wDiscriminatorDeep;
        const std::vector<double>&                      JetsAK8_hDiscriminatorDeep;
        const std::vector<int>&                         JetsAK8_multiplicity;
    
    JetAK8Collection(const NTupleReader& tr) 
            : JetsAK8(tr.getVec<TLorentzVector>("JetsAK8"))
            , JetsAK8_NsubjettinessTau1(tr.getVec<double>("JetsAK8_NsubjettinessTau1"))
            , JetsAK8_NsubjettinessTau2(tr.getVec<double>("JetsAK8_NsubjettinessTau2"))
            , JetsAK8_NsubjettinessTau3(tr.getVec<double>("JetsAK8_NsubjettinessTau3"))
            , JetsAK8_softDropMass(tr.getVec<double>("JetsAK8_softDropMass"))
            , JetsAK8_prunedMass(tr.getVec<double>("JetsAK8_prunedMass"))
            , JetsAK8_axismajor(tr.getVec<double>("JetsAK8_axismajor"))
            , JetsAK8_axisminor(tr.getVec<double>("JetsAK8_axisminor"))
            , JetsAK8_subjets(tr.getVec<std::vector<TLorentzVector>>("JetsAK8_subjets"))
            , JetsAK8_tDiscriminatorDeep(tr.getVec<double>("JetsAK8_tDiscriminatorDeep"))
            , JetsAK8_wDiscriminatorDeep(tr.getVec<double>("JetsAK8_wDiscriminatorDeep"))
            , JetsAK8_hDiscriminatorDeep(tr.getVec<double>("JetsAK8_hDiscriminatorDeep"))
            , JetsAK8_multiplicity(tr.getVec<int>("JetsAK8_multiplicity"))
        {
        }
    };

    class Factor
    {
    private:
        const std::vector<double>& Jets_jerFactor;
        const std::vector<double>& Jets_jecUnc;
        const std::vector<double>& JetsJECup_jerFactor;
        const std::vector<double>& JetsJECdown_jerFactor;
        const std::vector<double>& Jets_jerFactorUp;
        const std::vector<double>& Jets_jerFactorDown;

    public:
        double factor(const std::string& name, const int i, const int j) const
        {            
            if     (name == "JECup")   return (1./Jets_jerFactor.at(i))*(1+Jets_jecUnc.at(i))*JetsJECup_jerFactor.at(j);
            else if(name == "JECdown") return (1./Jets_jerFactor.at(i))*(1-Jets_jecUnc.at(i))*JetsJECdown_jerFactor.at(j);
            else if(name == "JERup")   return (1./Jets_jerFactor.at(i))*Jets_jerFactorUp.at(i);
            else if(name == "JERdown") return (1./Jets_jerFactor.at(i))*Jets_jerFactorDown.at(i);                
            else                       return -9999.9;
        }
        
        Factor(const NTupleReader& tr, const std::string& name = "") 
            : Jets_jerFactor(tr.getVec<double>("Jets"+name+"_jerFactor"))
            , Jets_jecUnc(tr.getVec<double>("Jets"+name+"_jecUnc"))
            , JetsJECup_jerFactor(tr.getVec<double>("Jets"+name+"JECup_jerFactor"))
            , JetsJECdown_jerFactor(tr.getVec<double>("Jets"+name+"JECdown_jerFactor"))
            , Jets_jerFactorUp(tr.getVec<double>("Jets"+name+"_jerFactorUp"))
            , Jets_jerFactorDown(tr.getVec<double>("Jets"+name+"_jerFactorDown"))
        {
        }
    };

    void deriveJetCollection(NTupleReader& tr, const JetCollection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("Jets"+name+"_origIndex");
               
        auto& newJets                             = tr.createDerivedVec<TLorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_bDiscriminatorCSV           = tr.createDerivedVec<double>("Jets"+name+"_bDiscriminatorCSV", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb         = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb        = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb          = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb    = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds   = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobc         = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobudsg      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobudsg", jc.Jets.size());
        auto& newJets_qgLikelihood                = tr.createDerivedVec<double>("Jets"+name+"_qgLikelihood", jc.Jets.size());
        auto& newJets_ID                          = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor                = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());
        auto& newJets_ptD                         = tr.createDerivedVec<double>("Jets"+name+"_ptD", jc.Jets.size());
        auto& newJets_axismajor                   = tr.createDerivedVec<double>("Jets"+name+"_axismajor", jc.Jets.size());
        auto& newJets_axisminor                   = tr.createDerivedVec<double>("Jets"+name+"_axisminor", jc.Jets.size());
        auto& newJets_multiplicity                = tr.createDerivedVec<int>("Jets"+name+"_multiplicity", jc.Jets.size());
        auto& newJets_neutralEmEnergyFraction     = tr.createDerivedVec<double>("Jets"+name+"_neutralEmEnergyFraction", jc.Jets.size());
        auto& newJets_chargedEmEnergyFraction     = tr.createDerivedVec<double>("Jets"+name+"_chargedEmEnergyFraction", jc.Jets.size());
        auto& newJets_neutralHadronEnergyFraction = tr.createDerivedVec<double>("Jets"+name+"_neutralHadronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronEnergyFraction = tr.createDerivedVec<double>("Jets"+name+"_chargedHadronEnergyFraction", jc.Jets.size());
        auto& newJets_muonEnergyFraction          = tr.createDerivedVec<double>("Jets"+name+"_muonEnergyFraction", jc.Jets.size());
        auto& newJets_hfHadronEnergyFraction      = tr.createDerivedVec<double>("Jets"+name+"_hfHadronEnergyFraction", jc.Jets.size());
        auto& newJets_hfEMEnergyFraction          = tr.createDerivedVec<double>("Jets"+name+"_hfEMEnergyFraction", jc.Jets.size());
        auto& newJets_photonEnergyFraction        = tr.createDerivedVec<double>("Jets"+name+"_photonEnergyFraction", jc.Jets.size());
        auto& newJets_electronEnergyFraction      = tr.createDerivedVec<double>("Jets"+name+"_electronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_chargedHadronMultiplicity", jc.Jets.size());
        auto& newJets_neutralHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_neutralHadronMultiplicity", jc.Jets.size());
        auto& newJets_photonMultiplicity          = tr.createDerivedVec<int>("Jets"+name+"_photonMultiplicity", jc.Jets.size());
        auto& newJets_electronMultiplicity        = tr.createDerivedVec<int>("Jets"+name+"_electronMultiplicity", jc.Jets.size());
        auto& newJets_muonMultiplicity            = tr.createDerivedVec<int>("Jets"+name+"_muonMultiplicity", jc.Jets.size());

        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i = newIndex[newJets_origIndex.at(j)];
            
            newJets.at(j)                             = jc.Jets.at(i)*f.factor(name,i,j);
            newJets_bDiscriminatorCSV.at(j)           = jc.Jets_bDiscriminatorCSV.at(i);
            newJets_bJetTagDeepCSVprobb.at(j)         = jc.Jets_bJetTagDeepCSVprobb.at(i);
            newJets_bJetTagDeepCSVprobbb.at(j)        = jc.Jets_bJetTagDeepCSVprobbb.at(i);
            newJets_bJetTagDeepCSVtotb.at(j)          = ( jc.Jets_bJetTagDeepCSVprobb.at(i) + jc.Jets_bJetTagDeepCSVprobbb.at(i) );
            newJets_bJetTagDeepFlavourprobb.at(j)     = jc.Jets_bJetTagDeepFlavourprobb.at(i);
            newJets_bJetTagDeepFlavourprobbb.at(j)    = jc.Jets_bJetTagDeepFlavourprobbb.at(i);
            newJets_bJetTagDeepFlavourtotb.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobb.at(i) + jc.Jets_bJetTagDeepFlavourprobbb.at(i) );
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
        }
    }

    void derivePtMassScaledJetCollection(NTupleReader& tr, const JetCollection& jc, const std::string& name, double scalePt, double scaleMass)
    {
        tr.registerDerivedVar("JetID"+name, jc.JetID);        

        auto& newJets                             = tr.createDerivedVec<TLorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_bDiscriminatorCSV           = tr.createDerivedVec<double>("Jets"+name+"_bDiscriminatorCSV", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb         = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb        = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb          = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb    = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc     = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds   = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());
        auto& newJets_qgLikelihood                = tr.createDerivedVec<double>("Jets"+name+"_qgLikelihood", jc.Jets.size());
        auto& newJets_ID                          = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor                = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());
        auto& newJets_ptD                         = tr.createDerivedVec<double>("Jets"+name+"_ptD", jc.Jets.size());
        auto& newJets_axismajor                   = tr.createDerivedVec<double>("Jets"+name+"_axismajor", jc.Jets.size());
        auto& newJets_axisminor                   = tr.createDerivedVec<double>("Jets"+name+"_axisminor", jc.Jets.size());
        auto& newJets_multiplicity                = tr.createDerivedVec<int>("Jets"+name+"_multiplicity", jc.Jets.size());
        auto& newJets_neutralEmEnergyFraction     = tr.createDerivedVec<double>("Jets"+name+"_neutralEmEnergyFraction", jc.Jets.size());
        auto& newJets_chargedEmEnergyFraction     = tr.createDerivedVec<double>("Jets"+name+"_chargedEmEnergyFraction", jc.Jets.size());
        auto& newJets_neutralHadronEnergyFraction = tr.createDerivedVec<double>("Jets"+name+"_neutralHadronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronEnergyFraction = tr.createDerivedVec<double>("Jets"+name+"_chargedHadronEnergyFraction", jc.Jets.size());
        auto& newJets_muonEnergyFraction          = tr.createDerivedVec<double>("Jets"+name+"_muonEnergyFraction", jc.Jets.size());
        auto& newJets_hfHadronEnergyFraction      = tr.createDerivedVec<double>("Jets"+name+"_hfHadronEnergyFraction", jc.Jets.size());
        auto& newJets_hfEMEnergyFraction          = tr.createDerivedVec<double>("Jets"+name+"_hfEMEnergyFraction", jc.Jets.size());
        auto& newJets_photonEnergyFraction        = tr.createDerivedVec<double>("Jets"+name+"_photonEnergyFraction", jc.Jets.size());
        auto& newJets_electronEnergyFraction      = tr.createDerivedVec<double>("Jets"+name+"_electronEnergyFraction", jc.Jets.size());
        auto& newJets_chargedHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_chargedHadronMultiplicity", jc.Jets.size());
        auto& newJets_neutralHadronMultiplicity   = tr.createDerivedVec<int>("Jets"+name+"_neutralHadronMultiplicity", jc.Jets.size());
        auto& newJets_photonMultiplicity          = tr.createDerivedVec<int>("Jets"+name+"_photonMultiplicity", jc.Jets.size());
        auto& newJets_electronMultiplicity        = tr.createDerivedVec<int>("Jets"+name+"_electronMultiplicity", jc.Jets.size());
        auto& newJets_muonMultiplicity            = tr.createDerivedVec<int>("Jets"+name+"_muonMultiplicity", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobc         = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobudsg      = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobudsg", jc.Jets.size());

        for(unsigned j = 0; j < jc.Jets.size(); ++j)
        {
            newJets[j].SetPtEtaPhiM( scalePt*jc.Jets[j].Pt(), jc.Jets[j].Eta(), jc.Jets[j].Phi(), scaleMass*jc.Jets[j].M() );
            newJets_bDiscriminatorCSV.at(j)           = jc.Jets_bDiscriminatorCSV.at(j);
            newJets_bJetTagDeepCSVprobb.at(j)         = jc.Jets_bJetTagDeepCSVprobb.at(j);
            newJets_bJetTagDeepCSVprobbb.at(j)        = jc.Jets_bJetTagDeepCSVprobbb.at(j);
            newJets_bJetTagDeepCSVtotb.at(j)          = ( jc.Jets_bJetTagDeepCSVprobb.at(j) + jc.Jets_bJetTagDeepCSVprobbb.at(j) );
            newJets_bJetTagDeepFlavourprobb.at(j)     = jc.Jets_bJetTagDeepFlavourprobb.at(j);
            newJets_bJetTagDeepFlavourprobbb.at(j)    = jc.Jets_bJetTagDeepFlavourprobbb.at(j);
            newJets_bJetTagDeepFlavourtotb.at(j)      = ( jc.Jets_bJetTagDeepFlavourprobb.at(j) + jc.Jets_bJetTagDeepFlavourprobbb.at(j) );
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
        }
    }

    void deriveJetAK8Collection(NTupleReader& tr, const JetAK8Collection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("JetsAK8"+name+"_origIndex");

        auto& newJetsAK8                    = tr.createDerivedVec<TLorentzVector>("JetsAK8"+name, jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau1  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau1", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau2  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau2", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau3  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau3", jc.JetsAK8.size());
        auto& newJetsAK8_softDropMass       = tr.createDerivedVec<double>("JetsAK8"+name+"_softDropMass", jc.JetsAK8.size());
        auto& newJetsAK8_prunedMass         = tr.createDerivedVec<double>("JetsAK8"+name+"_prunedMass", jc.JetsAK8.size());
        auto& newJetsAK8_axismajor          = tr.createDerivedVec<double>("JetsAK8"+name+"_axismajor", jc.JetsAK8.size());
        auto& newJetsAK8_axisminor          = tr.createDerivedVec<double>("JetsAK8"+name+"_axisminor", jc.JetsAK8.size());
        auto& newJetsAK8_subjets            = tr.createDerivedVec<std::vector<TLorentzVector>>("JetsAK8"+name+"_subjets", jc.JetsAK8.size());
        auto& newJetsAK8_tDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_tDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_wDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_wDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_hDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_hDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_multiplicity       = tr.createDerivedVec<int>("JetsAK8"+name+"_multiplicity", jc.JetsAK8.size());
             
        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i = newIndex[newJets_origIndex.at(j)];

            newJetsAK8.at(j)                    = jc.JetsAK8.at(i)*f.factor(name,i,j);
            newJetsAK8_NsubjettinessTau1.at(j)  = jc.JetsAK8_NsubjettinessTau1.at(i);
            newJetsAK8_NsubjettinessTau2.at(j)  = jc.JetsAK8_NsubjettinessTau2.at(i);
            newJetsAK8_NsubjettinessTau3.at(j)  = jc.JetsAK8_NsubjettinessTau3.at(i);
            newJetsAK8_softDropMass.at(j)       = jc.JetsAK8_softDropMass.at(i);
            newJetsAK8_prunedMass.at(j)         = jc.JetsAK8_prunedMass.at(i);
            newJetsAK8_axismajor.at(j)          = jc.JetsAK8_axismajor.at(i);
            newJetsAK8_axisminor.at(j)          = jc.JetsAK8_axisminor.at(i);
            newJetsAK8_subjets.at(j)            = jc.JetsAK8_subjets.at(i);
            newJetsAK8_tDiscriminatorDeep.at(j) = jc.JetsAK8_tDiscriminatorDeep.at(i);
            newJetsAK8_wDiscriminatorDeep.at(j) = jc.JetsAK8_wDiscriminatorDeep.at(i);
            newJetsAK8_hDiscriminatorDeep.at(j) = jc.JetsAK8_hDiscriminatorDeep.at(i);
            newJetsAK8_multiplicity.at(j)       = jc.JetsAK8_multiplicity.at(i);
        }
    }

    void derivePtMassScaledJetAK8Collection(NTupleReader& tr, const JetAK8Collection& jc, const std::string& name, double scalePt = 1.0, double scaleMass = 1.0)
    {
        auto& newJetsAK8                    = tr.createDerivedVec<TLorentzVector>("JetsAK8"+name, jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau1  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau1", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau2  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau2", jc.JetsAK8.size());
        auto& newJetsAK8_NsubjettinessTau3  = tr.createDerivedVec<double>("JetsAK8"+name+"_NsubjettinessTau3", jc.JetsAK8.size());
        auto& newJetsAK8_softDropMass       = tr.createDerivedVec<double>("JetsAK8"+name+"_softDropMass", jc.JetsAK8.size());
        auto& newJetsAK8_prunedMass         = tr.createDerivedVec<double>("JetsAK8"+name+"_prunedMass", jc.JetsAK8.size());
        auto& newJetsAK8_axismajor          = tr.createDerivedVec<double>("JetsAK8"+name+"_axismajor", jc.JetsAK8.size());
        auto& newJetsAK8_axisminor          = tr.createDerivedVec<double>("JetsAK8"+name+"_axisminor", jc.JetsAK8.size());
        auto& newJetsAK8_subjets            = tr.createDerivedVec<std::vector<TLorentzVector>>("JetsAK8"+name+"_subjets", jc.JetsAK8.size());
        auto& newJetsAK8_tDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_tDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_wDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_wDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_hDiscriminatorDeep = tr.createDerivedVec<double>("JetsAK8"+name+"_hDiscriminatorDeep", jc.JetsAK8.size());
        auto& newJetsAK8_multiplicity       = tr.createDerivedVec<int>("JetsAK8"+name+"_multiplicity", jc.JetsAK8.size());
             
        for(unsigned j = 0; j < jc.JetsAK8.size(); ++j)
        {
            newJetsAK8.at(j).SetPtEtaPhiM( scalePt*jc.JetsAK8[j].Pt(), jc.JetsAK8[j].Eta(), jc.JetsAK8[j].Phi(), scaleMass*jc.JetsAK8[j].M() );
            newJetsAK8_NsubjettinessTau1.at(j)  = jc.JetsAK8_NsubjettinessTau1.at(j);
            newJetsAK8_NsubjettinessTau2.at(j)  = jc.JetsAK8_NsubjettinessTau2.at(j);
            newJetsAK8_NsubjettinessTau3.at(j)  = jc.JetsAK8_NsubjettinessTau3.at(j);
            newJetsAK8_softDropMass.at(j)       = jc.JetsAK8_softDropMass.at(j);
            newJetsAK8_prunedMass.at(j)         = jc.JetsAK8_prunedMass.at(j);
            newJetsAK8_axismajor.at(j)          = jc.JetsAK8_axismajor.at(j);
            newJetsAK8_axisminor.at(j)          = jc.JetsAK8_axisminor.at(j);
            newJetsAK8_subjets.at(j)            = jc.JetsAK8_subjets.at(j);
            newJetsAK8_tDiscriminatorDeep.at(j) = jc.JetsAK8_tDiscriminatorDeep.at(j);
            newJetsAK8_wDiscriminatorDeep.at(j) = jc.JetsAK8_wDiscriminatorDeep.at(j);
            newJetsAK8_hDiscriminatorDeep.at(j) = jc.JetsAK8_hDiscriminatorDeep.at(j);
            newJetsAK8_multiplicity.at(j)       = jc.JetsAK8_multiplicity.at(j);
        }
    }

    void prepNTupleVars(NTupleReader& tr)
    {
        // Creating the jet pT and mass scaled collection
        JetCollection jc(tr);
        JetAK8Collection jcAK8(tr);
        const auto& runYear = tr.getVar<std::string>("runYear");
        derivePtMassScaledJetCollection(tr, jc, "mpTScaled", pTMass_[runYear].first, pTMass_[runYear].second);
        derivePtMassScaledJetAK8Collection(tr, jcAK8, "mpTScaled", 1.0, 1.0);

        // Create DeepCSV b-jet discriminator vector
        const auto& Jets_bJetTagDeepCSVprobb       = tr.getVec<double>("Jets_bJetTagDeepCSVprobb");
        const auto& Jets_bJetTagDeepCSVprobbb      = tr.getVec<double>("Jets_bJetTagDeepCSVprobbb");
              auto& Jets_bJetTagDeepCSVtotb        = tr.createDerivedVec<double>("Jets_bJetTagDeepCSVtotb", Jets_bJetTagDeepCSVprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobb   = tr.getVec<double>("Jets_bJetTagDeepFlavourprobb");
        const auto& Jets_bJetTagDeepFlavourprobbb  = tr.getVec<double>("Jets_bJetTagDeepFlavourprobbb");
              auto& Jets_bJetTagDeepFlavourtotb    = tr.createDerivedVec<double>("Jets_bJetTagDeepFlavourtotb", Jets_bJetTagDeepFlavourprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobc   = tr.getVec<double>("Jets_bJetTagDeepFlavourprobc");
        const auto& Jets_bJetTagDeepFlavourprobuds = tr.getVec<double>("Jets_bJetTagDeepFlavourprobuds");
              auto& Jets_bJetTagDeepFlavourtotq    = tr.createDerivedVec<double>("Jets_bJetTagDeepFlavourtotq", Jets_bJetTagDeepFlavourprobc.size());

        for(unsigned j = 0; j < Jets_bJetTagDeepCSVprobb.size(); ++j)
        {
            Jets_bJetTagDeepCSVtotb.at(j)     = ( Jets_bJetTagDeepCSVprobb.at(j) + Jets_bJetTagDeepCSVprobbb.at(j) );
            Jets_bJetTagDeepFlavourtotb.at(j) = ( Jets_bJetTagDeepFlavourprobb.at(j) + Jets_bJetTagDeepFlavourprobbb.at(j) );
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

        // Create the eventCounter variable to keep track of processed events
        int w = 1;
        if(runtype == "MC")
        {
            const auto& Weight = tr.getVar<double>("Weight");
            w = (Weight >= 0.0) ? 1 : -1;
        }
        tr.registerDerivedVar<int>("eventCounter",w);        
    }
public:
    PrepNTupleVars() : pTMass_({{"2016",{0.95,0.95}},{"2017",{0.95,1.01}},{"2018pre",{0.95,0.98}},{"2018post",{0.95,0.98}}})
    {
    }

    void operator()(NTupleReader& tr)
    {
        prepNTupleVars(tr);
    }
};

#endif
