#ifndef BTagCalibrationStandalone_H
#define BTagCalibrationStandalone_H

#include <TF1.h>
#include <TH1.h>
#include <map>
#include <vector>
#include <istream>
#include <ostream>
#include <memory>
#include <string>

class BTagEntry
{
public:
    enum OperatingPoint{ OP_LOOSE='L', OP_MEDIUM='M', OP_TIGHT='T', OP_RESHAPING='3' };
    enum JetFlavor{ FLAV_B=5, FLAV_C=4, FLAV_UDSG=0 };
    
    class Parameters 
    {
    public:
        OperatingPoint operatingPoint;
        std::string measurementType, sysType;
        JetFlavor jetFlavor;
        float etaMin, etaMax;
        float ptMin, ptMax;
        float discrMin, discrMax;
        
        Parameters(OperatingPoint op=OP_TIGHT, std::string measurement_type="comb", std::string sys_type="central", JetFlavor jf=FLAV_B,
                   float eta_min=-99999.0, float eta_max=99999.0, float pt_min=0.0, float pt_max=99999.0, float discr_min=0.0, float discr_max=99999.0);
    };
    
    BTagEntry() {}
    BTagEntry(const std::string& csvLine);
    BTagEntry(const std::string& func, Parameters p);
    BTagEntry(const TF1* func, Parameters p);
    BTagEntry(const TH1* histo, Parameters p);
    ~BTagEntry() {}
    static std::string makeCSVHeader();
    std::string makeCSVLine() const;
    static std::string trimStr(std::string str);
    std::string formula;
    Parameters params;    
};

/**
 * BTagCalibration
 *
 * The 'hierarchy' of stored information is this:
 * - by tagger (BTagCalibration)
 *   - by operating point or reshape bin
 *     - by jet parton flavor
 *       - by type of measurement
 *         - by systematic
 *           - by eta bin
 *             - as 1D-function dependent of pt or discriminant
 *
 ************************************************************/
class BTagCalibration
{
public:
    BTagCalibration() {}
    BTagCalibration(const std::string& tagger);
    BTagCalibration(const std::string& tagger, const std::string& filename);
    ~BTagCalibration() {}
    
    std::string tagger() const {return tagger_;}
    void addEntry(const BTagEntry& entry);
    const std::vector<BTagEntry>& getEntries(const BTagEntry::Parameters& par) const;
    void readCSV(std::istream& s);
    void readCSV(const std::string& s);
    void makeCSV(std::ostream& s) const;
    std::string makeCSV() const;
    
protected:
    static std::string token(const BTagEntry::Parameters& par);

    std::string tagger_;
    std::map<std::string, std::vector<BTagEntry>> data_;
};

/**
 * BTagCalibrationReader
 *
 * Helper class to pull out a specific set of BTagEntry's out of a
 * BTagCalibration. TF1 functions are set up at initialization time.
 *
 ************************************************************/
class BTagCalibrationReader
{
public:
    class BTagCalibrationReaderImpl
    {
        friend class BTagCalibrationReader;

    public:
        class TmpEntry 
        {
        public:
            float etaMin, etaMax, ptMin, ptMax, discrMin, discrMax;
            TF1 func;

            TmpEntry(float etaMin=-99999.0, float etaMax=99999.0, float ptMin=0.0, float ptMax=99999.0, float discrMin=0.0, float discrMax=99999.0, TF1 func=TF1());
        };

    private:
        BTagCalibrationReaderImpl(BTagEntry::OperatingPoint op, const std::string& sysType);
        BTagCalibrationReaderImpl(BTagEntry::OperatingPoint op, const std::string& sysType, const std::vector<std::string> & otherSysTypes);

        void load(const BTagCalibration& c, BTagEntry::JetFlavor jf, std::string measurementType);
        float eval(BTagEntry::JetFlavor jf, float eta, float pt, float discr) const;
        float eval_auto_bounds(const std::string& sys, BTagEntry::JetFlavor jf, float eta, float pt, float discr) const;

        std::pair<float, float> min_max_pt(BTagEntry::JetFlavor jf, float eta, float discr) const;
        BTagEntry::OperatingPoint op_;
        std::string sysType_;
        std::vector<std::vector<TmpEntry> > tmpData_;  // first index: jetFlavor
        std::vector<bool> useAbsEta_;                  // first index: jetFlavor
        std::map<std::string, BTagCalibrationReaderImpl*> otherSysTypeReaders_;
    };
    
    BTagCalibrationReader() {}
    BTagCalibrationReader(BTagEntry::OperatingPoint op, const std::string& sysType="central");
    BTagCalibrationReader(BTagEntry::OperatingPoint op, const std::string& sysType, const std::vector<std::string> & otherSysTypes);
    
    void load(const BTagCalibration& c, BTagEntry::JetFlavor jf, const std::string & measurementType="comb");
    float eval(BTagEntry::JetFlavor jf, float eta, float pt, float discr=0.0) const;
    float eval_auto_bounds(const std::string& sys, BTagEntry::JetFlavor jf, float eta, float pt, float discr=0.0) const;
    std::pair<float, float> min_max_pt(BTagEntry::JetFlavor jf, float eta, float discr=0.0) const;
    
protected:
    BTagCalibrationReaderImpl* pimpl;
};

#endif
