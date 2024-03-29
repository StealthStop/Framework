#include "../../../Framework/Framework/include/samples.h"

#include <string>
#include <iostream>
extern "C" {
    AnaSamples::SampleSet* SS_new(char *ssfile)
    {
        return new AnaSamples::SampleSet(ssfile);
    }
    AnaSamples::SampleCollection* SC_new(AnaSamples::SampleSet* ss, char* scfile)
    {
        return new AnaSamples::SampleCollection(scfile, *ss); 
    }
    AnaSamples::SampleCollection* SSSC_new(char *ssfile, char* scfile)
    {
        AnaSamples::SampleSet *ss = new AnaSamples::SampleSet(ssfile);
        return new AnaSamples::SampleCollection(scfile, *ss); 
    }
    int SC_samples_size(AnaSamples::SampleCollection* sc, char *scn){ return (*sc)[std::string(scn)].size(); }
    char const ** SC_samples(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = (*sc)[std::string(scn)];
        const char **array = new const char*[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            std::string* s = new std::string (sample.filePath + "/" + sample.fileName);
            array[i++] = s->c_str();
        }
        return array;
    }
    char const ** SC_samples_names(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = sc->getSampleLabels(std::string(scn));
        const char **array = new const char*[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            array[i++] = sample.c_str();
        }
        return array;
    }
    char const ** SS_samples_treePaths(AnaSamples::SampleSet* ss)
    {
        const char **array = new const char*[ss->size()];
        int i = 0;
        for(auto& sample : *ss)
        {
            array[i++] = sample.second.treePath.c_str();
        }
        return array;
    }

    int const * SC_samples_nGenEvts(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = (*sc)[std::string(scn)];
        int *array = new int[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            array[i++] = sample.nGenEvts;
        }
        return array;
    }
    int const * SC_samples_nActEvts(AnaSamples::SampleCollection* sc, char *scn)
    {
        auto& sampleVec = (*sc)[std::string(scn)];
        int *array = new int[sampleVec.size()];
        int i = 0;
        for(auto& sample : sampleVec)
        {
            array[i++] = sample.nActEvts;
        }
        return array;
    }

    char const ** SS_samples(AnaSamples::SampleSet* ss)
    {
        const char **array = new const char*[ss->size()];
        int i = 0;
        for(auto& sample : *ss)
        {
            std::string* s = new std::string (sample.second.filePath + "/" + sample.second.fileName);
            array[i++] = s->c_str();
        }
        return array;
    }
    char const ** SS_samples_names(AnaSamples::SampleSet* ss)
    {
        const char **array = new const char*[ss->size()];
        int i = 0;
        for(auto& sample : *ss)
        {
            array[i++] = sample.first.c_str();
        }
        return array;
    }
    int SC_samplecollection_size(AnaSamples::SampleCollection* sc, char *scn){ return sc->size(); }
    int SS_samples_size(AnaSamples::SampleCollection* ss){ return ss->size(); }
    char const ** SC_samplecollection_names(AnaSamples::SampleCollection* sc)
    {
        const char **array = new const char*[sc->size()];
        int i = 0;
        for(auto& sample : *sc)
        {
            array[i++] = sample.first.c_str();
        }
        return array;
    }
}
