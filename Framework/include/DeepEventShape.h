#ifndef DEEPEVENTSHAPE_H
#define DEEPEVENTSHAPE_H

#include "tensorflow/c/c_api.h"
#include "TopTagger/CfgParser/include/Context.hh"
#include "TopTagger/CfgParser/include/CfgDocument.hh"

#include "cstdlib"
#include "cstdio"
#include "cstring"

class EventShapeCalculator
{
private:
    float* basePtr_;
    std::string myVarSuffix_;
    int len_;

    int fwm2_top6_, fwm3_top6_, fwm4_top6_, fwm5_top6_, fwm6_top6_, fwm7_top6_, fwm8_top6_, fwm9_top6_, fwm10_top6_, jmt_ev0_top6_, jmt_ev1_top6_, jmt_ev2_top6_;
    int NGoodJets_double_;
    int Jet_pt_1_, Jet_pt_2_, Jet_pt_3_, Jet_pt_4_, Jet_pt_5_, Jet_pt_6_, Jet_pt_7_;
    int Jet_eta_1_, Jet_eta_2_, Jet_eta_3_, Jet_eta_4_, Jet_eta_5_, Jet_eta_6_, Jet_eta_7_;
    int Jet_phi_1_, Jet_phi_2_, Jet_phi_3_, Jet_phi_4_, Jet_phi_5_, Jet_phi_6_, Jet_phi_7_;
    int Jet_m_1_, Jet_m_2_, Jet_m_3_, Jet_m_4_, Jet_m_5_, Jet_m_6_, Jet_m_7_;
    int GoodLeptons_pt_1_, GoodLeptons_eta_1_, GoodLeptons_phi_1_, GoodLeptons_m_1_;
    int BestComboAvgMass_;

public:
    EventShapeCalculator(std::string myVarSuffix = "")
        : myVarSuffix_(myVarSuffix)
    {
        fwm2_top6_ = fwm3_top6_ = fwm4_top6_ = fwm5_top6_ = fwm6_top6_ = fwm7_top6_ = fwm8_top6_ = fwm9_top6_ = fwm10_top6_ = jmt_ev0_top6_ = jmt_ev1_top6_ = jmt_ev2_top6_ = -1; 
        NGoodJets_double_ = -1;
        Jet_pt_1_ = Jet_pt_2_ = Jet_pt_3_ = Jet_pt_4_ = Jet_pt_5_ = Jet_pt_6_ = Jet_pt_7_ = -1;
        Jet_eta_1_ = Jet_eta_2_ = Jet_eta_3_ = Jet_eta_4_ = Jet_eta_5_ = Jet_eta_6_ = Jet_eta_7_ = -1;
        Jet_phi_1_ = Jet_phi_2_ = Jet_phi_3_ = Jet_phi_4_ = Jet_phi_5_ = Jet_phi_6_ = Jet_phi_7_ = -1;
        Jet_m_1_ = Jet_m_2_ = Jet_m_3_ = Jet_m_4_ = Jet_m_5_ = Jet_m_6_ = Jet_m_7_ = -1;
        GoodLeptons_pt_1_ = GoodLeptons_eta_1_ = GoodLeptons_phi_1_ = GoodLeptons_m_1_ = -1;
        BestComboAvgMass_ = -1;
    }

    /**
     *The job of mapVars is to populate the internal offests for all variables in the input variable list with their memory location in the data array.  To be called only once.
     */
    void mapVars(const std::vector<std::string>& vars)
    {
        for(unsigned int i = 0; i < vars.size(); ++i)
        {
            if(     vars[i].compare("fwm2_top6") == 0)  fwm2_top6_ = i;
            else if(vars[i].compare("fwm3_top6") == 0)  fwm3_top6_ = i;
            else if(vars[i].compare("fwm4_top6") == 0)  fwm4_top6_ = i;
            else if(vars[i].compare("fwm5_top6") == 0)  fwm5_top6_ = i;
            else if(vars[i].compare("fwm6_top6") == 0)  fwm6_top6_ = i;
            else if(vars[i].compare("fwm7_top6") == 0)  fwm7_top6_ = i;
            else if(vars[i].compare("fwm8_top6") == 0)  fwm8_top6_ = i;
            else if(vars[i].compare("fwm9_top6") == 0)  fwm9_top6_ = i;
            else if(vars[i].compare("fwm10_top6") == 0) fwm10_top6_ = i;
            else if(vars[i].compare("jmt_ev0_top6") == 0) jmt_ev0_top6_ = i;
            else if(vars[i].compare("jmt_ev1_top6") == 0) jmt_ev1_top6_ = i;
            else if(vars[i].compare("jmt_ev2_top6") == 0) jmt_ev2_top6_ = i;
            else if(vars[i].compare("NGoodJets_double") == 0) NGoodJets_double_ = i;
            else if(vars[i].compare("Jet_pt_1") == 0) Jet_pt_1_ = i;
            else if(vars[i].compare("Jet_pt_2") == 0) Jet_pt_2_ = i;
            else if(vars[i].compare("Jet_pt_3") == 0) Jet_pt_3_ = i;
            else if(vars[i].compare("Jet_pt_4") == 0) Jet_pt_4_ = i;
            else if(vars[i].compare("Jet_pt_5") == 0) Jet_pt_5_ = i;
            else if(vars[i].compare("Jet_pt_6") == 0) Jet_pt_6_ = i;
            else if(vars[i].compare("Jet_pt_7") == 0) Jet_pt_7_ = i;
            else if(vars[i].compare("Jet_eta_1") == 0) Jet_eta_1_ = i;
            else if(vars[i].compare("Jet_eta_2") == 0) Jet_eta_2_ = i;
            else if(vars[i].compare("Jet_eta_3") == 0) Jet_eta_3_ = i;
            else if(vars[i].compare("Jet_eta_4") == 0) Jet_eta_4_ = i;
            else if(vars[i].compare("Jet_eta_5") == 0) Jet_eta_5_ = i;
            else if(vars[i].compare("Jet_eta_6") == 0) Jet_eta_6_ = i;
            else if(vars[i].compare("Jet_eta_7") == 0) Jet_eta_7_ = i;
            else if(vars[i].compare("Jet_phi_1") == 0) Jet_phi_1_ = i;
            else if(vars[i].compare("Jet_phi_2") == 0) Jet_phi_2_ = i;
            else if(vars[i].compare("Jet_phi_3") == 0) Jet_phi_3_ = i;
            else if(vars[i].compare("Jet_phi_4") == 0) Jet_phi_4_ = i;
            else if(vars[i].compare("Jet_phi_5") == 0) Jet_phi_5_ = i;
            else if(vars[i].compare("Jet_phi_6") == 0) Jet_phi_6_ = i;
            else if(vars[i].compare("Jet_phi_7") == 0) Jet_phi_7_ = i;
            else if(vars[i].compare("Jet_m_1") == 0) Jet_m_1_ = i;
            else if(vars[i].compare("Jet_m_2") == 0) Jet_m_2_ = i;
            else if(vars[i].compare("Jet_m_3") == 0) Jet_m_3_ = i;
            else if(vars[i].compare("Jet_m_4") == 0) Jet_m_4_ = i;
            else if(vars[i].compare("Jet_m_5") == 0) Jet_m_5_ = i;
            else if(vars[i].compare("Jet_m_6") == 0) Jet_m_6_ = i;
            else if(vars[i].compare("Jet_m_7") == 0) Jet_m_7_ = i;
            else if(vars[i].compare("GoodLeptons_pt_1") == 0) GoodLeptons_pt_1_ = i;
            else if(vars[i].compare("GoodLeptons_eta_1") == 0) GoodLeptons_eta_1_ = i;
            else if(vars[i].compare("GoodLeptons_phi_1") == 0) GoodLeptons_phi_1_ = i;
            else if(vars[i].compare("GoodLeptons_m_1") == 0) GoodLeptons_m_1_ = i;
            else if(vars[i].compare("BestComboAvgMass") == 0) BestComboAvgMass_ = i;
        }
    }
    /**
     *The job of setPtr is to set the starting place of memory block where the data will be written. To be called only once for the creation of the array pointed to by data.
     */
    void setPtr(float* data) {basePtr_ = data;}
    /**
     *Calculate the requested variables and store the values directly in the input array for the MVA
     */
    void calculateVars(const NTupleReader& tr)
    {
        if(fwm2_top6_ >= 0)  *(basePtr_ + fwm2_top6_) =  tr.getVar<double>("fwm2_top6"+myVarSuffix_);
        if(fwm3_top6_ >= 0)  *(basePtr_ + fwm3_top6_) =  tr.getVar<double>("fwm3_top6"+myVarSuffix_);
        if(fwm4_top6_ >= 0)  *(basePtr_ + fwm4_top6_) =  tr.getVar<double>("fwm4_top6"+myVarSuffix_);
        if(fwm5_top6_ >= 0)  *(basePtr_ + fwm5_top6_) =  tr.getVar<double>("fwm5_top6"+myVarSuffix_);
        if(fwm6_top6_ >= 0)  *(basePtr_ + fwm6_top6_) =  tr.getVar<double>("fwm6_top6"+myVarSuffix_);
        if(fwm7_top6_ >= 0)  *(basePtr_ + fwm7_top6_) =  tr.getVar<double>("fwm7_top6"+myVarSuffix_);
        if(fwm8_top6_ >= 0)  *(basePtr_ + fwm8_top6_) =  tr.getVar<double>("fwm8_top6"+myVarSuffix_);
        if(fwm9_top6_ >= 0)  *(basePtr_ + fwm9_top6_) =  tr.getVar<double>("fwm9_top6"+myVarSuffix_);
        if(fwm10_top6_ >= 0) *(basePtr_ + fwm10_top6_) =  tr.getVar<double>("fwm10_top6"+myVarSuffix_);
        if(jmt_ev0_top6_ >= 0) *(basePtr_ + jmt_ev0_top6_) =  tr.getVar<double>("jmt_ev0_top6"+myVarSuffix_);
        if(jmt_ev1_top6_ >= 0) *(basePtr_ + jmt_ev1_top6_) =  tr.getVar<double>("jmt_ev1_top6"+myVarSuffix_);
        if(jmt_ev2_top6_ >= 0) *(basePtr_ + jmt_ev2_top6_) =  tr.getVar<double>("jmt_ev2_top6"+myVarSuffix_);
        if(NGoodJets_double_ >= 0) *(basePtr_ + NGoodJets_double_) =  static_cast<double>(tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_));
        if(Jet_pt_1_  >= 0) *(basePtr_ + Jet_pt_1_ ) = tr.getVar<double>("Jet_pt_1"+myVarSuffix_);
        if(Jet_pt_2_  >= 0) *(basePtr_ + Jet_pt_2_ ) = tr.getVar<double>("Jet_pt_2"+myVarSuffix_);
        if(Jet_pt_3_  >= 0) *(basePtr_ + Jet_pt_3_ ) = tr.getVar<double>("Jet_pt_3"+myVarSuffix_);
        if(Jet_pt_4_  >= 0) *(basePtr_ + Jet_pt_4_ ) = tr.getVar<double>("Jet_pt_4"+myVarSuffix_);
        if(Jet_pt_5_  >= 0) *(basePtr_ + Jet_pt_5_ ) = tr.getVar<double>("Jet_pt_5"+myVarSuffix_);
        if(Jet_pt_6_  >= 0) *(basePtr_ + Jet_pt_6_ ) = tr.getVar<double>("Jet_pt_6"+myVarSuffix_);
        if(Jet_pt_7_  >= 0) *(basePtr_ + Jet_pt_7_ ) = tr.getVar<double>("Jet_pt_7"+myVarSuffix_);
        if(Jet_eta_1_  >= 0) *(basePtr_ + Jet_eta_1_ ) = tr.getVar<double>("Jet_eta_1"+myVarSuffix_);
        if(Jet_eta_2_  >= 0) *(basePtr_ + Jet_eta_2_ ) = tr.getVar<double>("Jet_eta_2"+myVarSuffix_);
        if(Jet_eta_3_  >= 0) *(basePtr_ + Jet_eta_3_ ) = tr.getVar<double>("Jet_eta_3"+myVarSuffix_);
        if(Jet_eta_4_  >= 0) *(basePtr_ + Jet_eta_4_ ) = tr.getVar<double>("Jet_eta_4"+myVarSuffix_);
        if(Jet_eta_5_  >= 0) *(basePtr_ + Jet_eta_5_ ) = tr.getVar<double>("Jet_eta_5"+myVarSuffix_);
        if(Jet_eta_6_  >= 0) *(basePtr_ + Jet_eta_6_ ) = tr.getVar<double>("Jet_eta_6"+myVarSuffix_);
        if(Jet_eta_7_  >= 0) *(basePtr_ + Jet_eta_7_ ) = tr.getVar<double>("Jet_eta_7"+myVarSuffix_);
        if(Jet_phi_1_  >= 0) *(basePtr_ + Jet_phi_1_ ) = tr.getVar<double>("Jet_phi_1"+myVarSuffix_);
        if(Jet_phi_2_  >= 0) *(basePtr_ + Jet_phi_2_ ) = tr.getVar<double>("Jet_phi_2"+myVarSuffix_);
        if(Jet_phi_3_  >= 0) *(basePtr_ + Jet_phi_3_ ) = tr.getVar<double>("Jet_phi_3"+myVarSuffix_);
        if(Jet_phi_4_  >= 0) *(basePtr_ + Jet_phi_4_ ) = tr.getVar<double>("Jet_phi_4"+myVarSuffix_);
        if(Jet_phi_5_  >= 0) *(basePtr_ + Jet_phi_5_ ) = tr.getVar<double>("Jet_phi_5"+myVarSuffix_);
        if(Jet_phi_6_  >= 0) *(basePtr_ + Jet_phi_6_ ) = tr.getVar<double>("Jet_phi_6"+myVarSuffix_);
        if(Jet_phi_7_  >= 0) *(basePtr_ + Jet_phi_7_ ) = tr.getVar<double>("Jet_phi_7"+myVarSuffix_);
        if(Jet_m_1_  >= 0) *(basePtr_ + Jet_m_1_ ) = tr.getVar<double>("Jet_m_1"+myVarSuffix_);
        if(Jet_m_2_  >= 0) *(basePtr_ + Jet_m_2_ ) = tr.getVar<double>("Jet_m_2"+myVarSuffix_);
        if(Jet_m_3_  >= 0) *(basePtr_ + Jet_m_3_ ) = tr.getVar<double>("Jet_m_3"+myVarSuffix_);
        if(Jet_m_4_  >= 0) *(basePtr_ + Jet_m_4_ ) = tr.getVar<double>("Jet_m_4"+myVarSuffix_);
        if(Jet_m_5_  >= 0) *(basePtr_ + Jet_m_5_ ) = tr.getVar<double>("Jet_m_5"+myVarSuffix_);
        if(Jet_m_6_  >= 0) *(basePtr_ + Jet_m_6_ ) = tr.getVar<double>("Jet_m_6"+myVarSuffix_);
        if(Jet_m_7_  >= 0) *(basePtr_ + Jet_m_7_ ) = tr.getVar<double>("Jet_m_7"+myVarSuffix_);
        if(GoodLeptons_pt_1_  >= 0) *(basePtr_ + GoodLeptons_pt_1_ ) = tr.getVar<double>("GoodLeptons_pt_1"+myVarSuffix_);
        if(GoodLeptons_eta_1_ >= 0) *(basePtr_ + GoodLeptons_eta_1_) = tr.getVar<double>("GoodLeptons_eta_1"+myVarSuffix_);
        if(GoodLeptons_phi_1_ >= 0) *(basePtr_ + GoodLeptons_phi_1_) = tr.getVar<double>("GoodLeptons_phi_1"+myVarSuffix_);
        if(GoodLeptons_m_1_   >= 0) *(basePtr_ + GoodLeptons_m_1_  ) = tr.getVar<double>("GoodLeptons_m_1"+myVarSuffix_);
        if(BestComboAvgMass_ >= 0) *(basePtr_ + BestComboAvgMass_) = tr.getVar<double>("BestComboAvgMass"+myVarSuffix_);
    }
};

class DeepEventShape
{
private:
    double discriminator_;
    std::string modelFile_, inputOp_, outputOp_, myVarSuffix_;
    int minNJet_, maxNJet_;

    //Tensoflow session pointer
    TF_Session* session_;

    //Input variable names 
    std::vector<std::string> vars_;
    std::vector<double> binEdges_;

    std::vector<TF_Output>     inputs_;
    std::vector<TF_Output>     outputs_;
    std::vector<TF_Operation*> targets_;

    //variable calclator
    std::shared_ptr<EventShapeCalculator> varCalculator_;

    template<typename T> std::vector<T> getVecFromCfg(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& var, const cfg::Context& localCxt, const T& defaultYo)
    {
        std::vector<T> vec;
        int iVar = 0;
        bool keepLooping;
        do
        {
            keepLooping = false;
        
            //Get variable name
            T v = cfgDoc->get(var, iVar, localCxt, defaultYo);
        
            //if it is a non empty string save in vector
            if(v != defaultYo)
            {
                keepLooping = true;
        
                vec.push_back(v);
            }
            ++iVar;
        }
        while(keepLooping);
        
        return vec;
    }

    void getParameters(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& localContextName)
    {
        //Construct contexts
        cfg::Context localCxt(localContextName);

        modelFile_   = cfgDoc->get("modelFile", localCxt, "");
        inputOp_     = cfgDoc->get("inputOp",   localCxt, "main_input");
        outputOp_    = cfgDoc->get("outputOp",  localCxt, "first_output/Softmax");
        minNJet_     = cfgDoc->get("minNJet",   localCxt, 7);
        maxNJet_     = cfgDoc->get("maxNJet",   localCxt, 7);
        vars_        = getVecFromCfg<std::string>(cfgDoc, "mvaVar", localCxt, "");
        binEdges_    = getVecFromCfg<double>(cfgDoc, "binEdges", localCxt, -1);
        
        //Variable to hold tensorflow status
        TF_Status* status = TF_NewStatus();
        
        //get the grafdef from the file
        TF_Buffer* graph_def = read_file(modelFile_);
        
        // Import graph_def into graph
        TF_Graph* graph = TF_NewGraph();
        TF_ImportGraphDefOptions* graph_opts = TF_NewImportGraphDefOptions();
        TF_GraphImportGraphDef(graph, graph_def, graph_opts, status);
        TF_DeleteImportGraphDefOptions(graph_opts);
        TF_DeleteBuffer(graph_def);
        
        //Create tensorflow session from imported graph
        TF_SessionOptions* sess_opts = TF_NewSessionOptions();
        uint8_t config[] = {0x10, 0x01};
        TF_SetConfig(sess_opts, static_cast<void*>(config), 2, status);
        session_ = TF_NewSession(graph, sess_opts, status);
        TF_DeleteSessionOptions(sess_opts);
        
        TF_Operation* op_x = TF_GraphOperationByName(graph, inputOp_.c_str());
        TF_Operation* op_y = TF_GraphOperationByName(graph, outputOp_.c_str());
        
        //Clean up graph
        TF_DeleteGraph(graph);
        
        inputs_ .emplace_back(TF_Output({op_x, 0}));
        outputs_.emplace_back(TF_Output({op_y, 0}));
        targets_.emplace_back(op_y);
        
        TF_DeleteStatus(status);

        //map variables
        varCalculator_.reset(new EventShapeCalculator(myVarSuffix_));
        varCalculator_->mapVars(vars_);
    }

    void runDeepEventShape(NTupleReader& tr)
    {
        //tensorflow status variable
        TF_Status* status = TF_NewStatus();
        
        //Create place to store the output vectors 
        std::vector<TF_Tensor*>    output_values(1);
        
        //Construct tensorflow input tensor
        std::vector<TF_Tensor*> input_values;
        const int elemSize = sizeof(float);
        std::vector<int64_t> dims = {static_cast<int64_t>(1), static_cast<int64_t>(vars_.size())};
        int nelem = 1;
        for(const auto dimLen : dims) nelem *= dimLen;
        TF_Tensor* input_values_0 =  TF_AllocateTensor(TF_FLOAT, dims.data(), dims.size(), elemSize*nelem);
        
        input_values = { input_values_0 };
        varCalculator_->setPtr(static_cast<float*>(TF_TensorData(input_values_0)));

        varCalculator_->calculateVars(tr);

        //predict values
        TF_SessionRun(session_,
                      // RunOptions
                      nullptr,
                      // Input tensors
                      inputs_.data(), input_values.data(), inputs_.size(),
                      // Output tensors
                      outputs_.data(), output_values.data(), outputs_.size(),
                      // Target operations
                      targets_.data(), targets_.size(),
                      // RunMetadata
                      nullptr,
                      // Output status
                      status);
        
        //Get output discriminators 
        auto discriminators = static_cast<float*>(TF_TensorData(output_values[0]));                
        
        //discriminators is a 2D array, we only want the first entry of every array
        double discriminator = static_cast<double>(discriminators[0]);

        for(auto* tensor : input_values)  TF_DeleteTensor(tensor);
        for(auto* tensor : output_values) TF_DeleteTensor(tensor);
        
        TF_DeleteStatus(status);

        // Register Variables
        tr.registerDerivedVar("deepESM_val"+myVarSuffix_, discriminator);

        // Define and register deepESM bins
        const auto& NGoodJets_pt30 = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        int nMVABin = (binEdges_.size() / (maxNJet_ - minNJet_ + 1)) - 1;
        int nJetBinning;
        if(NGoodJets_pt30 < minNJet_) nJetBinning = 0;
        else if(minNJet_ <= NGoodJets_pt30 && NGoodJets_pt30 <= maxNJet_) nJetBinning = NGoodJets_pt30-minNJet_;
        else if(maxNJet_ < NGoodJets_pt30) nJetBinning = maxNJet_-minNJet_;

        for(int i = (nMVABin+1)*nJetBinning + 1; i < (nMVABin+1)*(nJetBinning+1); i++)
        {
            bool passDeepESMBin = discriminator > binEdges_[i-1] && discriminator <= binEdges_[i];
            int bin = i - (nMVABin+1)*nJetBinning;
            tr.registerDerivedVar("deepESM_bin"+std::to_string(bin)+myVarSuffix_, passDeepESMBin);
            //std::cout<<"nMVABin: "<<nMVABin<<" NJets: "<<NGoodJets_pt30<<" nJetBinning: "<<nJetBinning
            //         <<" i: "<<i<<" lowBinEdge: "<<binEdges_[i-1]<<" highBinEdge: "<<binEdges_[i]<<" MVABinNumber: "<<bin<<std::endl;
        }
    }

    static void free_buffer(void* data, size_t length) 
    {
        free(data);
    }

    TF_Buffer* read_file(const std::string& file) 
    {
        FILE* f = fopen(file.c_str(), "rb");

        fseek(f, 0, SEEK_END);
        long fsize = ftell(f);
        fseek(f, 0, SEEK_SET);  //same as rewind(f);

        void* data = malloc(fsize);
        fread(data, fsize, 1, f);
        fclose(f);

        TF_Buffer* buf = TF_NewBuffer();
        buf->data = data;
        buf->length = fsize;
        buf->data_deallocator = free_buffer;
        return buf;
    }

public:
    DeepEventShape(DeepEventShape&& husk) 
        : discriminator_(husk.discriminator_)
        , modelFile_(husk.modelFile_)
        , inputOp_(husk.inputOp_)
        , outputOp_(husk.outputOp_)
        , myVarSuffix_(husk.myVarSuffix_)
        , minNJet_(husk.minNJet_)
        , maxNJet_(husk.maxNJet_)
        , session_(husk.session_)
        , vars_(husk.vars_)
        , binEdges_(husk.binEdges_)
        , inputs_(husk.inputs_)
        , outputs_(husk.outputs_)
        , targets_(husk.targets_)
        , varCalculator_(husk.varCalculator_)
    {
        husk.session_ = nullptr;
    }
    
    DeepEventShape(const std::string cfgFileName = "DeepEventShape.cfg", std::string localContextName = "Info", bool printStatus = true, const std::string myVarSuffix = "")
        : myVarSuffix_(myVarSuffix)
    {
        if(printStatus) std::cout<<"Setting up DeepEventShape"<<std::endl;
        
        //buffer to hold file contents 
        std::string cfgText;

        FILE *f = fopen(cfgFileName.c_str(), "r");
        char buff[1024];
        for(; !feof(f) && fgets(buff, 1023, f);)
        {
            cfgText += buff;
        }
        
        fclose(f);
        
        //pass raw text to cfg parser, to return parsed document
        std::unique_ptr<cfg::CfgDocument> cfgDoc = cfg::CfgDocument::parseDocument(cfgText);
        getParameters(cfgDoc, localContextName);

        if(printStatus) std::cout<<"Using "+cfgFileName+" as the DeepEventShape config file"<<std::endl;
    }

    ~DeepEventShape()
    {
        if(session_)
        {
            TF_Status* status = TF_NewStatus();
            TF_CloseSession(session_, status);
            TF_DeleteSession(session_, status);
            TF_DeleteStatus(status);
        }
    }
    
    void operator()(NTupleReader& tr)
    {
        runDeepEventShape(tr);
    }
};

#endif
