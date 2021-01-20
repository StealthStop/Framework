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
        const std::vector<double>& Jets_bJetTagDeepCSVprobb;
        const std::vector<double>& Jets_bJetTagDeepCSVprobbb;
        const std::vector<double>& Jets_bJetTagDeepFlavourprobb;
        const std::vector<double>& Jets_bJetTagDeepFlavourprobbb;
        const std::vector<double>& Jets_bJetTagDeepFlavourprobc;
        const std::vector<double>& Jets_bJetTagDeepFlavourprobuds;
        const std::vector<double>& Jets_bJetTagDeepFlavourprobg;
        const std::vector<bool>& Jets_ID;
        const bool& JetID;
        const std::vector<int>& Jets_partonFlavor;
        
        JetCollection(const NTupleReader& tr) 
            : Jets(tr.getVec<TLorentzVector>("Jets"))
            , Jets_bJetTagDeepCSVprobb(tr.getVec<double>("Jets_bJetTagDeepCSVprobb"))
            , Jets_bJetTagDeepCSVprobbb(tr.getVec<double>("Jets_bJetTagDeepCSVprobbb"))
            , Jets_bJetTagDeepFlavourprobb(tr.getVec<double>("Jets_bJetTagDeepFlavourprobb"))
            , Jets_bJetTagDeepFlavourprobbb(tr.getVec<double>("Jets_bJetTagDeepFlavourprobbb"))
            , Jets_bJetTagDeepFlavourprobc(tr.getVec<double>("Jets_bJetTagDeepFlavourprobc"))
            , Jets_bJetTagDeepFlavourprobuds(tr.getVec<double>("Jets_bJetTagDeepFlavourprobuds"))
            , Jets_bJetTagDeepFlavourprobg(tr.getVec<double>("Jets_bJetTagDeepFlavourprobg"))
            , Jets_ID(tr.getVec<bool>("Jets_ID"))
            , JetID(tr.getVar<bool>("JetID"))
            , Jets_partonFlavor(tr.getVec<int>("Jets_partonFlavor"))
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
        
        Factor(const NTupleReader& tr) 
            : Jets_jerFactor(tr.getVec<double>("Jets_jerFactor"))
            , Jets_jecUnc(tr.getVec<double>("Jets_jecUnc"))
            , JetsJECup_jerFactor(tr.getVec<double>("JetsJECup_jerFactor"))
            , JetsJECdown_jerFactor(tr.getVec<double>("JetsJECdown_jerFactor"))
            , Jets_jerFactorUp(tr.getVec<double>("Jets_jerFactorUp"))
            , Jets_jerFactorDown(tr.getVec<double>("Jets_jerFactorDown"))
        {
        }
    };

    void deriveJetCollection(NTupleReader& tr, const JetCollection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("Jets"+name+"_origIndex");
               
        auto& newJets = tr.createDerivedVec<TLorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());

        auto& newJets_ID = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());
        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i = newIndex[newJets_origIndex.at(j)];
            
            newJets.at(j) = jc.Jets.at(i)*f.factor(name,i,j);
            newJets_bJetTagDeepCSVprobb.at(j) = jc.Jets_bJetTagDeepCSVprobb.at(i);
            newJets_bJetTagDeepCSVprobbb.at(j) = jc.Jets_bJetTagDeepCSVprobbb.at(i);
            newJets_bJetTagDeepCSVtotb.at(j) = ( jc.Jets_bJetTagDeepCSVprobb.at(i) + jc.Jets_bJetTagDeepCSVprobbb.at(i) );
            newJets_bJetTagDeepFlavourprobb.at(j) = jc.Jets_bJetTagDeepFlavourprobb.at(i);
            newJets_bJetTagDeepFlavourprobbb.at(j) = jc.Jets_bJetTagDeepFlavourprobbb.at(i);
            newJets_bJetTagDeepFlavourtotb.at(j)  = ( jc.Jets_bJetTagDeepFlavourprobb.at(i) + jc.Jets_bJetTagDeepFlavourprobbb.at(i) );
            newJets_bJetTagDeepFlavourtotq.at(j) = ( jc.Jets_bJetTagDeepFlavourprobc.at(i) + jc.Jets_bJetTagDeepFlavourprobuds.at(i) );
            newJets_bJetTagDeepFlavourprobg.at(j) = jc.Jets_bJetTagDeepFlavourprobg.at(i); 
            newJets_bJetTagDeepFlavourprobc.at(j) = jc.Jets_bJetTagDeepFlavourprobc.at(i); 
            newJets_bJetTagDeepFlavourprobuds.at(j) = jc.Jets_bJetTagDeepFlavourprobuds.at(i); 

            newJets_ID.at(j) = jc.Jets_ID.at(i);
            newJets_partonFlavor.at(j) = jc.Jets_partonFlavor.at(i);
        }
    }

    void derivePtMassScaledJetCollection(NTupleReader& tr, const JetCollection& jc, const std::string& name, double scalePt, double scaleMass)
    {
        tr.registerDerivedVar("JetID"+name, jc.JetID);        

        auto& newJets = tr.createDerivedVec<TLorentzVector>("Jets"+name, jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobbb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotb", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourtotq = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourtotq", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobg = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobg", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobc = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobc", jc.Jets.size());
        auto& newJets_bJetTagDeepFlavourprobuds = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepFlavourprobuds", jc.Jets.size());

        auto& newJets_ID = tr.createDerivedVec<bool>("Jets"+name+"_ID", jc.Jets.size());
        auto& newJets_partonFlavor = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor", jc.Jets.size());

        for(unsigned j = 0; j < jc.Jets.size(); ++j)
        {
            newJets[j].SetPtEtaPhiM( scalePt*jc.Jets[j].Pt(), jc.Jets[j].Eta(), jc.Jets[j].Phi(), scaleMass*jc.Jets[j].M() );
            newJets_bJetTagDeepCSVprobb.at(j) = jc.Jets_bJetTagDeepCSVprobb.at(j);
            newJets_bJetTagDeepCSVprobbb.at(j) = jc.Jets_bJetTagDeepCSVprobbb.at(j);
            newJets_bJetTagDeepCSVtotb.at(j) = ( jc.Jets_bJetTagDeepCSVprobb.at(j) + jc.Jets_bJetTagDeepCSVprobbb.at(j) );
            newJets_bJetTagDeepFlavourprobb.at(j) = jc.Jets_bJetTagDeepFlavourprobb.at(j);
            newJets_bJetTagDeepFlavourprobbb.at(j) = jc.Jets_bJetTagDeepFlavourprobbb.at(j);
            newJets_bJetTagDeepFlavourtotb.at(j)  = ( jc.Jets_bJetTagDeepFlavourprobb.at(j) + jc.Jets_bJetTagDeepFlavourprobbb.at(j) );
            newJets_bJetTagDeepFlavourtotq.at(j) = ( jc.Jets_bJetTagDeepFlavourprobc.at(j) + jc.Jets_bJetTagDeepFlavourprobuds.at(j) );
            newJets_bJetTagDeepFlavourprobg.at(j) = jc.Jets_bJetTagDeepFlavourprobg.at(j); 
            newJets_bJetTagDeepFlavourprobc.at(j) = jc.Jets_bJetTagDeepFlavourprobc.at(j); 
            newJets_bJetTagDeepFlavourprobuds.at(j) = jc.Jets_bJetTagDeepFlavourprobuds.at(j); 

            newJets_ID.at(j) = jc.Jets_ID.at(j);
            newJets_partonFlavor.at(j) = jc.Jets_partonFlavor.at(j);
        }
    }

    void prepNTupleVars(NTupleReader& tr)
    {
        // Creating the jet pT and mass scaled collection
        JetCollection jc(tr);
        const auto& runYear = tr.getVar<std::string>("runYear");
        derivePtMassScaledJetCollection(tr, jc, "mpTScaled", pTMass_[runYear].first, pTMass_[runYear].second);

        // Create DeepCSV b-jet discriminator vector
        const auto& Jets_bJetTagDeepCSVprobb  = tr.getVec<double>("Jets_bJetTagDeepCSVprobb");
        const auto& Jets_bJetTagDeepCSVprobbb = tr.getVec<double>("Jets_bJetTagDeepCSVprobbb");
        auto& Jets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets_bJetTagDeepCSVtotb", Jets_bJetTagDeepCSVprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobb  = tr.getVec<double>("Jets_bJetTagDeepFlavourprobb");
        const auto& Jets_bJetTagDeepFlavourprobbb = tr.getVec<double>("Jets_bJetTagDeepFlavourprobbb");
        auto& Jets_bJetTagDeepFlavourtotb = tr.createDerivedVec<double>("Jets_bJetTagDeepFlavourtotb", Jets_bJetTagDeepFlavourprobbb.size());
        const auto& Jets_bJetTagDeepFlavourprobc  = tr.getVec<double>("Jets_bJetTagDeepFlavourprobc");
        const auto& Jets_bJetTagDeepFlavourprobuds = tr.getVec<double>("Jets_bJetTagDeepFlavourprobuds");
        auto& Jets_bJetTagDeepFlavourtotq = tr.createDerivedVec<double>("Jets_bJetTagDeepFlavourtotq", Jets_bJetTagDeepFlavourprobc.size());
        for(unsigned j = 0; j < Jets_bJetTagDeepCSVprobb.size(); ++j)
        {
            Jets_bJetTagDeepCSVtotb.at(j) = ( Jets_bJetTagDeepCSVprobb.at(j) + Jets_bJetTagDeepCSVprobbb.at(j) );
            Jets_bJetTagDeepFlavourtotb.at(j)  = ( Jets_bJetTagDeepFlavourprobb.at(j) + Jets_bJetTagDeepFlavourprobbb.at(j) );
            Jets_bJetTagDeepFlavourtotq.at(j) = ( Jets_bJetTagDeepFlavourprobc.at(j) + Jets_bJetTagDeepFlavourprobuds.at(j) );
        }
 
        // Create JEC/R variation variables if needed
        const auto& runtype = tr.getVar<std::string>("runtype");
        if( !tr.checkBranchInTree("JetsJECup") && runtype == "MC")
        {
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
