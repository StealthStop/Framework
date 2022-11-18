#ifndef MinTupleMaker_h
#define MinTupleMaker_h

#include "NTupleReader/include/NTupleReader.h"

#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <set>
#include <map>
#include <string>

#include <iostream>

class MiniTupleMaker
{
public:
    MiniTupleMaker(TTree *);

    MiniTupleMaker(const std::string&, const std::string& = "tree");

    ~MiniTupleMaker();

    void setTupleVars(const std::set<std::string>&);

    // To use derived variables initBranches must be called after the first tuple event is read
    void initBranches(const NTupleReader&);

    void fill();
    void fill(const NTupleReader&);

private:
    TFile* const file_;
    TTree* const tree_;

    bool startedFilling_ = false;

    std::set<std::string> tupleVars_;

    // These collections will appear to be vector<float> when MiniTupleMaker checks their type
    // However, the type is actually ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> 
    // So this exception is here to correctly set the type when assigning a pointer
    std::set<std::string> LVexceptions_ = {"Jets",         "JetsAK8",  "JetsAK8_subjets", 
                                           "Electrons",    "Muons",    "Photons", 
                                           "GenElectrons", "GenMuons", "GenTaus", "GenParticles"};

    template<typename T> void prepVar(const NTupleReader& tr, const std::string& name)
    {
        TBranch *tb = tree_->GetBranch(name.c_str());
        if(!tb) tree_->Branch(name.c_str(), static_cast<T*>(const_cast<void*>(tr.getPtr<T>(name))));
        else       tb->SetAddress(const_cast<void*>(tr.getPtr<T>(name)));
    }

    template<typename T> void prepVec(const NTupleReader& tr, const std::string& name)
    {
        TBranch *tb = tree_->GetBranch(name.c_str());

        if(!tb) tree_->Branch(name.c_str(), static_cast<std::vector<T>**>(const_cast<void*>(tr.getVecPtr<T>(name))));
        else       tb->SetAddress(const_cast<void*>(tr.getVecPtr<T>(name)));
    }
};

#endif
