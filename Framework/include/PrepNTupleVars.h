#ifndef PREPNTUPLEVARS_H
#define PREPNTUPLEVARS_H

class PrepNTupleVars
{
private:
    class JetCollection
    {
    public:
        const std::vector<TLorentzVector>& Jets;
        const std::vector<double>& Jets_bDiscriminatorCSV;
        const std::vector<bool>& Jets_ID;
        const std::vector<int>& Jets_partonFlavor;
        
        JetCollection(const NTupleReader& tr) 
            : Jets(tr.getVec<TLorentzVector>("Jets"))
            , Jets_bDiscriminatorCSV(tr.getVec<double>("Jets_bDiscriminatorCSV"))
            , Jets_ID(tr.getVec<bool>("Jets_ID"))
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
        const double factor(const std::string& name, const int i, const int j) const
        {            
            if     (name == "JECup")   return (1./Jets_jerFactor.at(i))*(1+Jets_jecUnc.at(i))*JetsJECup_jerFactor.at(j);
            else if(name == "JECdown") return (1./Jets_jerFactor.at(i))*(1-Jets_jecUnc.at(i))*JetsJECdown_jerFactor.at(j);
            else if(name == "JERup")   return (1./Jets_jerFactor.at(i))*Jets_jerFactorUp.at(i);
            else if(name == "JERdown") return (1./Jets_jerFactor.at(i))*Jets_jerFactorDown.at(i);                
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
            
        auto& newJets = tr.createDerivedVec<TLorentzVector>("Jets"+name);
        auto& newJets_bDiscriminatorCSV = tr.createDerivedVec<double>("Jets"+name+"_bDiscriminatorCSV");
        auto& newJets_ID = tr.createDerivedVec<bool>("Jets"+name+"_ID");
        auto& newJets_partonFlavor = tr.createDerivedVec<int>("Jets"+name+"_partonFlavor");
        for(unsigned j = 0; j < newJets_origIndex.size(); ++j)
        {
            int i = newIndex[newJets_origIndex.at(j)];
            
            newJets.emplace_back( jc.Jets.at(i)*f.factor(name,i,j) );
            newJets_bDiscriminatorCSV.emplace_back( jc.Jets_bDiscriminatorCSV.at(i) );
            newJets_ID.emplace_back( jc.Jets_ID.at(i) );
            newJets_partonFlavor.emplace_back( jc.Jets_partonFlavor.at(i) );
        }
    }

    void prepNTupleVars(NTupleReader& tr)
    {
        const auto& runYear = tr.getVar<std::string>("runYear");
        
        if(runYear == "2017")
        {
            const auto& Jets_origIndex = tr.getVec<int>("Jets_origIndex");
            std::vector<int> newIndex(Jets_origIndex.size());
            for(unsigned j = 0; j < Jets_origIndex.size(); ++j)
            {
                newIndex[Jets_origIndex.at(j)] = j;
            }

            JetCollection jc(tr);
            Factor f(tr);
            deriveJetCollection(tr, jc, f, newIndex, "JECup");
            deriveJetCollection(tr, jc, f, newIndex, "JECdown");
            deriveJetCollection(tr, jc, f, newIndex, "JERup");
            deriveJetCollection(tr, jc, f, newIndex, "JERdown");
        }
    }
public:
    PrepNTupleVars()
    {
        std::cout<<"Preparing variables from the NTuples"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        prepNTupleVars(tr);
    }
};

#endif
