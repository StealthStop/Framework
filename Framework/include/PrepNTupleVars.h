#ifndef PREPNTUPLEVARS_H
#define PREPNTUPLEVARS_H

#include "Framework/Framework/include/Utility.h"

class PrepNTupleVars
{
private:
    bool printJetpTScaled;

    class JetCollection
    {
    public:
        const std::vector<TLorentzVector>& Jets;
        const std::vector<double>& Jets_bJetTagDeepCSVprobb;
        const std::vector<double>& Jets_bJetTagDeepCSVprobbb;
        const std::vector<bool>& Jets_ID;
        const bool& JetID;
        const std::vector<int>& Jets_partonFlavor;
        
        JetCollection(const NTupleReader& tr) 
            : Jets(tr.getVec<TLorentzVector>("Jets"))
            , Jets_bJetTagDeepCSVprobb(tr.getVec<double>("Jets_bJetTagDeepCSVprobb"))
            , Jets_bJetTagDeepCSVprobbb(tr.getVec<double>("Jets_bJetTagDeepCSVprobbb"))
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

    template<typename T> T* scaleJetPt(const T& Jets, double scale)
    {
        auto* jets = new T(Jets.size());
        for(unsigned int i = 0; i < Jets.size(); i++)
        {                
            (*jets)[i].SetPtEtaPhiM( scale*Jets[i].Pt(), Jets[i].Eta(), Jets[i].Phi(), Jets[i].M() );
        }
        return jets;
    }

    void deriveJetCollection(NTupleReader& tr, const JetCollection& jc, const Factor& f, const std::vector<int>& newIndex, const std::string& name)
    {
        const auto& newJets_origIndex = tr.getVec<int>("Jets"+name+"_origIndex");
            
        auto& newJets = tr.createDerivedVec<TLorentzVector>("Jets"+name);
        auto& newJets_bJetTagDeepCSVprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb");
        auto& newJets_bJetTagDeepCSVprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb");
        auto& newJets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb");
        auto& newJets_ID = tr.createDerivedVec<bool>("Jets"+name+"_ID");
        auto& newJets_partonFlavor = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor");
        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i = newIndex[newJets_origIndex.at(j)];
            
            newJets.emplace_back( jc.Jets.at(i)*f.factor(name,i,j) );
            newJets_bJetTagDeepCSVprobb.emplace_back( jc.Jets_bJetTagDeepCSVprobb.at(i) );
            newJets_bJetTagDeepCSVprobbb.emplace_back( jc.Jets_bJetTagDeepCSVprobbb.at(i) );
            newJets_bJetTagDeepCSVtotb.emplace_back( jc.Jets_bJetTagDeepCSVprobb.at(i) + jc.Jets_bJetTagDeepCSVprobbb.at(i) );
            newJets_ID.emplace_back( jc.Jets_ID.at(i) );
            newJets_partonFlavor.emplace_back( jc.Jets_partonFlavor.at(i) );
        }
    }

    void derivePtScaledJetCollection(NTupleReader& tr, const JetCollection& jc, const std::string& name, double scale)
    {
        auto* jets = scaleJetPt(jc.Jets, scale);
        tr.registerDerivedVec("Jets"+name,jets);
                
        auto& newJets_bJetTagDeepCSVprobb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobb");
        auto& newJets_bJetTagDeepCSVprobbb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVprobbb");
        auto& newJets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets"+name+"_bJetTagDeepCSVtotb");
        auto& newJets_ID = tr.createDerivedVec<bool>("Jets"+name+"_ID");
        tr.registerDerivedVar("JetID"+name, jc.JetID);
        auto& newJets_partonFlavor = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor");
        
        for(unsigned j = 0; j < jc.Jets.size(); ++j)
        {
            newJets_bJetTagDeepCSVprobb.emplace_back( jc.Jets_bJetTagDeepCSVprobb.at(j) );
            newJets_bJetTagDeepCSVprobbb.emplace_back( jc.Jets_bJetTagDeepCSVprobbb.at(j) );
            newJets_bJetTagDeepCSVtotb.emplace_back( jc.Jets_bJetTagDeepCSVprobb.at(j) + jc.Jets_bJetTagDeepCSVprobbb.at(j) );
            newJets_ID.emplace_back( jc.Jets_ID.at(j) );
            newJets_partonFlavor.emplace_back( jc.Jets_partonFlavor.at(j) );
        }
    }

    void prepNTupleVars(NTupleReader& tr)
    {
        // Reweighting the nominal jet collection (needs to be done first if doing it this way)
        double scale = 0.95;
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& scaleJetpT = tr.getVar<bool>("scaleJetpT");
        if(runtype == "MC" && scaleJetpT)
        {
            printJetpTScaled = true;
            auto& Jets = tr.getVec<TLorentzVector>("Jets");
            auto* jets = scaleJetPt(Jets, scale);
            tr.registerDerivedVec("Jets",jets);
        }

        // Creating the jet pT scaled collection
        JetCollection jc(tr);
        derivePtScaledJetCollection(tr, jc, "pTScaled", scale);

        // Create DeepCSV b-jet discriminator vector
        const auto& Jets_bJetTagDeepCSVprobb  = tr.getVec<double>("Jets_bJetTagDeepCSVprobb");
        const auto& Jets_bJetTagDeepCSVprobbb = tr.getVec<double>("Jets_bJetTagDeepCSVprobbb");
        auto& Jets_bJetTagDeepCSVtotb = tr.createDerivedVec<double>("Jets_bJetTagDeepCSVtotb");
        for(unsigned j = 0; j < Jets_bJetTagDeepCSVprobb.size(); ++j)
        {
            Jets_bJetTagDeepCSVtotb.emplace_back( Jets_bJetTagDeepCSVprobb.at(j) + Jets_bJetTagDeepCSVprobbb.at(j) );
        }
 
        // Create JEC/R variation variables if needed
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
    PrepNTupleVars() : printJetpTScaled(false)
    {
        std::cout<<"Preparing variables from the NTuples"<<std::endl;
    }

    ~PrepNTupleVars()
    {
        if(printJetpTScaled) std::cout<<utility::color("Warning: Scaled all Jet pT by 0.95","red")<<std::endl;
    }
    void operator()(NTupleReader& tr)
    {
        prepNTupleVars(tr);
    }
};

#endif
