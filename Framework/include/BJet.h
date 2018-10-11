#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::string myVarSuffix_;

    void setVar(bool pass, std::vector<bool>* boolVec, int& n)
    {
        if( pass )
        {
            boolVec->push_back(true);                
            n++;
        }
        else
        {
            boolVec->push_back(false);
        }            
    }

    void bjet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets"+myVarSuffix_+"_bDiscriminatorCSV");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& JetsID = tr.getVec<bool>("Jets"+myVarSuffix_+"_ID");
        const auto& GoodJets = tr.getVec<bool>("GoodJets");

        auto* bjets_loose_ = new std::vector<bool>();
        auto* bjets_pt30_loose_ = new std::vector<bool>();
        int NBJets_loose = 0, NBJets_pt30_loose = 0;

        auto* bjets_ = new std::vector<bool>();
        auto* bjets_pt30_ = new std::vector<bool>();
        auto* bjets_pt40_ = new std::vector<bool>();
        auto* bjets_pt45_ = new std::vector<bool>();
        int NBJets = 0, NBJets_pt30 = 0, NBJets_pt40 = 0, NBJets_pt45 = 0;

        auto* bjets_tight_ = new std::vector<bool>();
        auto* bjets_pt30_tight_ = new std::vector<bool>();
        auto* bjets_pt45_tight_ = new std::vector<bool>();
        int NBJets_tight = 0, NBJets_pt30_tight = 0, NBJets_pt45_tight = 0;


        auto* goodbjets_loose_ = new std::vector<bool>();
        auto* goodbjets_pt30_loose_ = new std::vector<bool>();
        int NGoodBJets_loose = 0, NGoodBJets_pt30_loose = 0;

        auto* goodbjets_ = new std::vector<bool>();
        auto* goodbjets_pt30_ = new std::vector<bool>();
        auto* goodbjets_pt40_ = new std::vector<bool>();
        auto* goodbjets_pt45_ = new std::vector<bool>();
        int NGoodBJets = 0, NGoodBJets_pt30 = 0, NGoodBJets_pt40 = 0, NGoodBJets_pt45 = 0;

        auto* goodbjets_tight_ = new std::vector<bool>();
        auto* goodbjets_pt30_tight_ = new std::vector<bool>();
        auto* goodbjets_pt45_tight_ = new std::vector<bool>();
        int NGoodBJets_tight = 0, NGoodBJets_pt30_tight = 0, NGoodBJets_pt45_tight = 0;

        double loose = 0.5426;
        double medium = 0.8484;
        double tight = 0.9535;

        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            TLorentzVector lv = Jets.at(ijet);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > loose                , bjets_loose_,      NBJets_loose     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > loose && lv.Pt() > 30, bjets_pt30_loose_, NBJets_pt30_loose);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > loose                 && GoodJets.at(ijet), goodbjets_loose_,      NGoodBJets_loose     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > loose && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_loose_, NGoodBJets_pt30_loose);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium                , bjets_,      NBJets     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 30, bjets_pt30_, NBJets_pt30);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 40, bjets_pt40_, NBJets_pt40);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 45, bjets_pt45_, NBJets_pt45);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium                 && GoodJets.at(ijet), goodbjets_,      NGoodBJets     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_, NGoodBJets_pt30);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 40 && GoodJets.at(ijet), goodbjets_pt40_, NGoodBJets_pt40);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > medium && lv.Pt() > 45 && GoodJets.at(ijet), goodbjets_pt45_, NGoodBJets_pt45);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight                , bjets_tight_,      NBJets_tight     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight && lv.Pt() > 30, bjets_pt30_tight_, NBJets_pt30_tight);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight && lv.Pt() > 45, bjets_pt45_tight_, NBJets_pt45_tight);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight                 && GoodJets.at(ijet), goodbjets_tight_,      NGoodBJets_tight     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_tight_, NGoodBJets_pt30_tight);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && Jets_bDiscriminatorCSV.at(ijet) > tight && lv.Pt() > 45 && GoodJets.at(ijet), goodbjets_pt45_tight_, NGoodBJets_pt45_tight);            
        }

        tr.registerDerivedVec("BJets_loose",        bjets_loose_);
        tr.registerDerivedVar("NBJets_loose",       NBJets_loose);
        tr.registerDerivedVec("BJets_pt30_loose",   bjets_pt30_loose_);
        tr.registerDerivedVar("NBJets_pt30_loose",  NBJets_pt30_loose);

        tr.registerDerivedVec("BJets",        bjets_);
        tr.registerDerivedVar("NBJets",       NBJets);
        tr.registerDerivedVec("BJets_pt30",   bjets_pt30_);
        tr.registerDerivedVar("NBJets_pt30",  NBJets_pt30);
        tr.registerDerivedVec("BJets_pt40",   bjets_pt40_);
        tr.registerDerivedVar("NBJets_pt40",  NBJets_pt40);
        tr.registerDerivedVec("BJets_pt45",   bjets_pt45_);
        tr.registerDerivedVar("NBJets_pt45",  NBJets_pt45);

        tr.registerDerivedVec("BJets_tight",        bjets_tight_);
        tr.registerDerivedVar("NBJets_tight",       NBJets_tight);
        tr.registerDerivedVec("BJets_pt30_tight",   bjets_pt30_tight_);
        tr.registerDerivedVar("NBJets_pt30_tight",  NBJets_pt30_tight);
        tr.registerDerivedVec("BJets_pt45_tight",   bjets_pt45_tight_);
        tr.registerDerivedVar("NBJets_pt45_tight",  NBJets_pt45_tight);

        tr.registerDerivedVec("GoodBJets_loose",        goodbjets_loose_);
        tr.registerDerivedVar("NGoodBJets_loose",       NGoodBJets_loose);
        tr.registerDerivedVec("GoodBJets_pt30_loose",   goodbjets_pt30_loose_);
        tr.registerDerivedVar("NGoodBJets_pt30_loose",  NGoodBJets_pt30_loose);

        tr.registerDerivedVec("GoodBJets",        goodbjets_);
        tr.registerDerivedVar("NGoodBJets",       NGoodBJets);
        tr.registerDerivedVec("GoodBJets_pt30",   goodbjets_pt30_);
        tr.registerDerivedVar("NGoodBJets_pt30",  NGoodBJets_pt30);
        tr.registerDerivedVec("GoodBJets_pt40",   goodbjets_pt40_);
        tr.registerDerivedVar("NGoodBJets_pt40",  NGoodBJets_pt40);
        tr.registerDerivedVec("GoodBJets_pt45",   goodbjets_pt45_);
        tr.registerDerivedVar("NGoodBJets_pt45",  NGoodBJets_pt45);
        
        tr.registerDerivedVec("GoodBJets_tight",        goodbjets_tight_);
        tr.registerDerivedVar("NGoodBJets_tight",       NGoodBJets_tight);
        tr.registerDerivedVec("GoodBJets_pt30_tight",   goodbjets_pt30_tight_);
        tr.registerDerivedVar("NGoodBJets_pt30_tight",  NGoodBJets_pt30_tight);
        tr.registerDerivedVec("GoodBJets_pt45_tight",   goodbjets_pt45_tight_);
        tr.registerDerivedVar("NGoodBJets_pt45_tight",  NGoodBJets_pt45_tight);
    }

public:
    BJet( std::string myVarSuffix = "" ) 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up BJet"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
