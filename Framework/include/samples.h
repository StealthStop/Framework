#ifndef ANASAMPLES_SAMPLES_H
#define ANASAMPLES_SAMPLES_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>

#include <cstring>
#include <algorithm>

namespace AnaSamples
{
  enum COLORS{kRed = 632, kGreen = 416, kBlack = 1, kMagenta = 616, kBlue = 600, kYellow = 400, kTeal = 840, kPink = 900, kOrange = 800, kSpring = 820, kWhite = 0, kGray = 0, kCyan = 432, kAzure = 860, kViolet = 880};

  class FileSummary
  {
   public:
    std::string tag;
    std::string filePath, fileName, treePath;
    double xsec, kfactor, nGenEvts, nActEvts;
    int color;
    bool isData_;
        
    FileSummary() {}
    FileSummary(const std::string& tag, const std::string& filePath, const std::string& fileName, const std::string& treePath, double xsec, double nGenEvts, double nActEvts, double kfactor, int color = kBlack) : tag(tag), filePath(filePath), fileName(fileName), treePath(treePath), xsec(xsec), kfactor(kfactor), nGenEvts(nGenEvts), nActEvts(nActEvts), color(color), isData_(false)
    {
      weight_ = xsec * kfactor / nGenEvts;
    }

    FileSummary(const std::string& tag, const std::string& filePath, const std::string& fileName, const std::string& treePath, double xsec, double nEvts, double kfactor, int color = kBlack) : tag(tag), filePath(filePath), fileName(fileName), treePath(treePath), xsec(xsec), kfactor(kfactor), nGenEvts(nEvts), nActEvts(nEvts), color(color), isData_(false)
    {
      weight_ = xsec * kfactor / nEvts;
    }

    //Constructor which doesn't make a xsec*kfactor/nEvts weighted sample, e.g. for use with data.
    //Initialize xsec, nEvts to 1 so that the comparison operators still work
    FileSummary(const std::string& tag, const std::string& filePath, const std::string& fileName, const std::string& treePath, double kfactor, int color = kBlack) : tag(tag), filePath(filePath), fileName(fileName), treePath(treePath), xsec(1), kfactor(kfactor), nGenEvts(1), nActEvts(1), color(color), isData_(true)
    {
      weight_ = kfactor;
    }

    double getWeight() const {return weight_;}

    const std::vector<std::string>& getFilelist() const {return filelist_;}
    template<class T> void addFilesToChain(T* chain, int startfile=0, int filerun=-1) const
    {
      if(filelist_.size() == 0) readFileList();
      if(filerun<0)filerun=filelist_.size();
      for(unsigned int fn = startfile; static_cast<int>(fn)<startfile+filerun && fn<filelist_.size(); fn++)
      {
        //printf("fn = %d, filelist_[fn]=%s\n", fn, filelist_[fn].c_str()); // testing
        chain->Add(filelist_[fn].c_str());
      }
    }
    mutable std::vector<std::string> filelist_;
    void addCollection(const std::string&);
    const std::set<std::string>& getCollections() const
    {
      return collections_;
    }
    void readFileList() const;

   private:
    double weight_;
    std::set<std::string> collections_;
  };

  bool operator< (const FileSummary& lhs, const FileSummary& rhs);
  bool operator== (const FileSummary& lhs, const FileSummary& rhs);
  bool operator!= (const FileSummary& lhs, const FileSummary& rhs);

  template<class T>
  class SampleBase
  {
   public:
    SampleBase() : nullT_() {};
    const T& operator[](const std::string& key) const
    {
            auto iter = sampleSet_.find(key);
            if(iter != sampleSet_.end()) return iter->second;
            else                         return nullT_;
    }
    const T& null() const {return nullT_;}
    
   protected:
    std::map<std::string, T> sampleSet_;
    const T nullT_;
    static constexpr unsigned int BUF_LEN_ = 4096;

    virtual bool parseCfgLine(const char* buf) = 0;

    void readCfg(const std::string& file)
    {
        FILE *fin = fopen(file.c_str(), "r");

        if(fin != nullptr)
        {
            char buf[BUF_LEN_];
            int lineNum = 0;
            while(!feof(fin) && fgets(buf, BUF_LEN_-1, fin))
            {
                ++lineNum;

                //skip comments
                if(strlen(buf) <= 0 || buf[0] == '#') continue;

                //skip empty lines (is there a way not to copy?)
                std::string sbuf(buf);
                if(std::all_of(sbuf.cbegin(), sbuf.cend(), [](char c) { return std::isspace(c); })) continue;

                //strip out commas
                char *bufPtr = strchr(buf, ',');
                while(bufPtr)
                {
                    *bufPtr = ' ';
                    bufPtr = strchr(bufPtr + 1, ',');
                }

                //strip out newlines
                bufPtr = strchr(buf, '\n');
                if(bufPtr) *bufPtr = '\0';

                //parse line
                if(!parseCfgLine(buf))
                {
                    //if false is returned, the line was not parsed properly
                    printf("Malformed line: %s:%i\n", file.c_str(), lineNum);
                }
            }

            fclose(fin);
        }
        else
        {
            std::cout << "ERROR: Unable to open file " << file << std::endl;
        }
    }


   public:
    decltype(sampleSet_.cbegin()) begin() const { return sampleSet_.cbegin(); }
    decltype(sampleSet_.cend())     end() const { return sampleSet_.cend(); }
    decltype(sampleSet_.size())    size() const { return sampleSet_.size(); }
  };

  class SampleCollection;
  
  class SampleSet : public SampleBase<FileSummary>
  {
    friend class SampleCollection;
   
   public:
    SampleSet(std::string file = "sampleSets.cfg", bool isCondor = false);
    void addSample(const std::string& tag, const std::string& filePath, const std::string& fileName, const std::string& treePath, double xsec, double nGenEvts, double nActEvts, double kfactor, int color = kBlack) 
    {
        sampleSet_[tag] = FileSummary(tag, filePath, fileName, treePath, xsec, nGenEvts, nActEvts, kfactor, color);
    }

    void addSample(const std::string& tag, const std::string& filePath, const std::string& fileName, const std::string& treePath, double xsec, double nEvts, double kfactor, int color = kBlack) 
    {
        sampleSet_[tag] = FileSummary(tag, filePath, fileName, treePath, xsec, nEvts, kfactor, color);
    }

    void addSample(const std::string& tag, const std::string& filePath, const std::string& fileName,  const std::string& treePath, double kfactor, int color = kBlack) 
    {
        sampleSet_[tag] = FileSummary(tag, filePath, fileName, treePath, kfactor, color);
    }
    
    // modify weights to compare two MC samples
    double getCrossSectionRatio(const std::vector<std::string>& sampleTags1, const std::vector<std::string>& sampleTags2, bool verbose = false);
    

   private:
    std::string fDir_;
    bool isCondor_;
    
    std::map<std::string, FileSummary>& getMap() { return sampleSet_; }
    
    bool parseCfgLine(const char* buf);
  };

  class SampleCollection : public SampleBase<std::vector<FileSummary>>
  {
   public:
    SampleCollection(const std::string& file, SampleSet& samples);
    std::vector<std::string>& getSampleLabels(std::string name);
    
    // modify weights to compare two MC samples
    double getCrossSectionRatio(std::string& sampleTag1, std::string sampleTag2, bool verbose = false);

   private:
    std::map<std::string, std::vector<std::string>> nameVec_;
    SampleSet& ss_;
    void addSampleSet(SampleSet& samples, const std::string& name, const std::vector<std::string>& vss);
    bool parseCfgLine(const char* buf);
  };
}

#endif
