#include "ExploreBackground.h"
#include "ExploreTopTagger.h"
#include "ExploreEventSelection.h"
#include "samples.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include<iostream>
#include <getopt.h>

int main(int argc, char *argv[])
{

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"doBackground",       no_argument, 0, 'b'},
        {"doTopTagger",        no_argument, 0, 't'},
        {"doEventSelection",   no_argument, 0, 's'},
        {"condor",           no_argument, 0, 'c'},
        {"histFile",   required_argument, 0, 'H'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
    };

    bool doBackground = false, doTopTagger = false, doEventSelection = false;
    bool runOnCondor = false;
    std::string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir;
    int nFiles = -1, startFile = 0, maxEvts = -1;

    while((opt = getopt_long(argc, argv, "btscH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'b':
            doBackground = true;
            break;

        case 't':
            doTopTagger = true;
            break;

        case 's':
            doEventSelection = true;
            break;

        case 'c':
            runOnCondor = true;
            break;

        case 'H':
            histFile = optarg;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            maxEvts = int(atoi(optarg));
            break;

        }
    }


    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        sampleloc = "condor";
    }

    AnaSamples::SampleSet        ss(sampleloc);
    AnaSamples::SampleCollection sc(ss);

    std::map<std::string, std::vector<AnaSamples::FileSummary>> fileMap;
    if(ss[dataSets] != ss.null())
    {
        fileMap[dataSets] = {ss[dataSets]};
        for(const auto& colls : ss[dataSets].getCollections())
        {
            fileMap[colls] = {ss[dataSets]};
        }
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = {sc[dataSets]};
        int i = 0;
        for(const auto& fs : sc[dataSets])
        {
            fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }
    std::set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    TFile* myfile = TFile::Open(histFile.c_str(), "RECREATE");

    TChain* ch = new TChain( (AnaSamples::treeName).c_str() ) ;

    if(doBackground)
    {
        ExploreBackground t = ExploreBackground(ch);
        std::cout << "Initializing..." << std::endl;
        t.InitHistos();
        for(const AnaSamples::FileSummary& file : vvf)
        {
            std::cout << "Running over sample " << file.tag << std::endl;
            TChain* new_ch = new TChain( (AnaSamples::treeName).c_str());
            t.Init(new_ch);
            file.addFilesToChain(new_ch, startFile, nFiles);
            double weight = file.getWeight();
            std::cout << "starting loop" << std::endl;
            std::string runtype = "";
            if(file.tag.find("Data") != std::string::npos)
                runtype = "Data";
            t.Loop(weight, maxEvts, runtype, file.tag);            
        }
        std::cout << "Writing histograms..." << std::endl;
        t.WriteHistos();
    }

    if(doTopTagger)
    {
        ExploreTopTagger t = ExploreTopTagger(ch);
        std::cout << "Initializing..." << std::endl;
        t.InitHistos();
        for(const AnaSamples::FileSummary& file : vvf)
        {
            std::cout << "Running over sample " << file.tag << std::endl;
            TChain* new_ch = new TChain( (AnaSamples::treeName).c_str());
            t.Init(new_ch);
            file.addFilesToChain(new_ch, startFile, nFiles);
            double weight = file.getWeight();
            std::string type = "";
            if(file.tag.find("qcd") != std::string::npos)
                type = "qcd";
            std::cout << "starting loop" << std::endl;
            t.Loop(type, weight, maxEvts);            
        }
        std::cout << "Writing histograms..." << std::endl;
        t.WriteHistos();
    }

    if(doEventSelection)
    {
        ExploreEventSelection t = ExploreEventSelection(ch);
        std::cout << "Initializing..." << std::endl;
        t.InitHistos();
        for(const AnaSamples::FileSummary& file : vvf)
        {
            std::cout << "Running over sample " << file.tag << std::endl;
            TChain* new_ch = new TChain( (AnaSamples::treeName).c_str());
            t.Init(new_ch);
            file.addFilesToChain(new_ch, startFile, nFiles);
            double weight = file.getWeight();
            std::cout << "starting loop" << std::endl;
            std::string runtype = "";
            if(file.tag.find("Data") != std::string::npos)
                runtype = "Data";
            t.Loop(weight, maxEvts, runtype, file.tag);
        }
        std::cout << "Writing histograms..." << std::endl;
        t.WriteHistos();
    }

    myfile->Close();

    return 0;
}
