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
    int len_;

    int fwm2_top6_, fwm3_top6_, fwm4_top6_, fwm5_top6_, fwm6_top6_, fwm7_top6_, fwm8_top6_, fwm9_top6_, fwm10_top6_, jmt_ev0_top6_, jmt_ev1_top6_, jmt_ev2_top6_;
    int NGoodJets_double_;
    int Jet_pt_1_, Jet_pt_2_, Jet_pt_3_, Jet_pt_4_, Jet_pt_5_, Jet_pt_6_, Jet_eta_1_, Jet_eta_2_, Jet_eta_3_, Jet_eta_4_, Jet_eta_5_, Jet_eta_6_;
    int Jet_phi_1_, Jet_phi_2_, Jet_phi_3_, Jet_phi_4_, Jet_phi_5_, Jet_phi_6_, Jet_m_1_, Jet_m_2_, Jet_m_3_, Jet_m_4_, Jet_m_5_, Jet_m_6_; 
    int BestComboAvgMass_;

public:
    EventShapeCalculator()
    {
        fwm2_top6_ = fwm3_top6_ = fwm4_top6_ = fwm5_top6_ = fwm6_top6_ = fwm7_top6_ = fwm8_top6_ = fwm9_top6_ = fwm10_top6_ = jmt_ev0_top6_ = jmt_ev1_top6_ = jmt_ev2_top6_ = -1; 
        NGoodJets_double_ = -1;
        Jet_pt_1_ = Jet_pt_2_ = Jet_pt_3_ = Jet_pt_4_ = Jet_pt_5_ = Jet_pt_6_ = Jet_eta_1_ = Jet_eta_2_ = Jet_eta_3_ = Jet_eta_4_ = Jet_eta_5_ = Jet_eta_6_ = -1;
        Jet_phi_1_ = Jet_phi_2_ = Jet_phi_3_ = Jet_phi_4_ = Jet_phi_5_ = Jet_phi_6_ = Jet_m_1_ = Jet_m_2_ = Jet_m_3_ = Jet_m_4_ = Jet_m_5_ = Jet_m_6_ = -1; 
        BestComboAvgMass_ = -1;
    }

    /**
     *The job of mapVars is to populate the internal offests for all variables in the input variable list with their memory location in the data array.  To be called only once.
     */
    void mapVars(const std::vector<std::string>& vars)
    {
        len_ = vars.size();

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
            else if(vars[i].compare("Jet_eta_1") == 0) Jet_eta_1_ = i;
            else if(vars[i].compare("Jet_eta_2") == 0) Jet_eta_2_ = i;
            else if(vars[i].compare("Jet_eta_3") == 0) Jet_eta_3_ = i;
            else if(vars[i].compare("Jet_eta_4") == 0) Jet_eta_4_ = i;
            else if(vars[i].compare("Jet_eta_5") == 0) Jet_eta_5_ = i;
            else if(vars[i].compare("Jet_eta_6") == 0) Jet_eta_6_ = i;
            else if(vars[i].compare("Jet_phi_1") == 0) Jet_phi_1_ = i;
            else if(vars[i].compare("Jet_phi_2") == 0) Jet_phi_2_ = i;
            else if(vars[i].compare("Jet_phi_3") == 0) Jet_phi_3_ = i;
            else if(vars[i].compare("Jet_phi_4") == 0) Jet_phi_4_ = i;
            else if(vars[i].compare("Jet_phi_5") == 0) Jet_phi_5_ = i;
            else if(vars[i].compare("Jet_phi_6") == 0) Jet_phi_6_ = i;
            else if(vars[i].compare("Jet_m_1") == 0) Jet_m_1_ = i;
            else if(vars[i].compare("Jet_m_2") == 0) Jet_m_2_ = i;
            else if(vars[i].compare("Jet_m_3") == 0) Jet_m_3_ = i;
            else if(vars[i].compare("Jet_m_4") == 0) Jet_m_4_ = i;
            else if(vars[i].compare("Jet_m_5") == 0) Jet_m_5_ = i;
            else if(vars[i].compare("Jet_m_6") == 0) Jet_m_6_ = i;
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
    void calculateVars(const NTupleReader& tr, int iCand)
    {
        if(fwm2_top6_ >= 0)  *(basePtr_ + fwm2_top6_ + len_*iCand) =  tr.getVar<double>("fwm2_top6");
        if(fwm3_top6_ >= 0)  *(basePtr_ + fwm3_top6_ + len_*iCand) =  tr.getVar<double>("fwm3_top6");
        if(fwm4_top6_ >= 0)  *(basePtr_ + fwm4_top6_ + len_*iCand) =  tr.getVar<double>("fwm4_top6");
        if(fwm5_top6_ >= 0)  *(basePtr_ + fwm5_top6_ + len_*iCand) =  tr.getVar<double>("fwm5_top6");
        if(fwm6_top6_ >= 0)  *(basePtr_ + fwm6_top6_ + len_*iCand) =  tr.getVar<double>("fwm6_top6");
        if(fwm7_top6_ >= 0)  *(basePtr_ + fwm7_top6_ + len_*iCand) =  tr.getVar<double>("fwm7_top6");
        if(fwm8_top6_ >= 0)  *(basePtr_ + fwm8_top6_ + len_*iCand) =  tr.getVar<double>("fwm8_top6");
        if(fwm9_top6_ >= 0)  *(basePtr_ + fwm9_top6_ + len_*iCand) =  tr.getVar<double>("fwm9_top6");
        if(fwm10_top6_ >= 0) *(basePtr_ + fwm10_top6_ + len_*iCand) =  tr.getVar<double>("fwm10_top6");
        if(jmt_ev0_top6_ >= 0) *(basePtr_ + jmt_ev0_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev0_top6");
        if(jmt_ev1_top6_ >= 0) *(basePtr_ + jmt_ev1_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev1_top6");
        if(jmt_ev2_top6_ >= 0) *(basePtr_ + jmt_ev2_top6_ + len_*iCand) =  tr.getVar<double>("jmt_ev2_top6");
        if(NGoodJets_double_ >= 0) *(basePtr_ + NGoodJets_double_ + len_*iCand) =  static_cast<double>(tr.getVar<unsigned long>("NGoodJets"));
        if(Jet_pt_1_ >= 0) *(basePtr_ + Jet_pt_1_ + len_*iCand) = tr.getVar<double>("Jet_pt_1");
        if(Jet_pt_2_ >= 0) *(basePtr_ + Jet_pt_2_ + len_*iCand) = tr.getVar<double>("Jet_pt_2");
        if(Jet_pt_3_ >= 0) *(basePtr_ + Jet_pt_3_ + len_*iCand) = tr.getVar<double>("Jet_pt_3");
        if(Jet_pt_4_ >= 0) *(basePtr_ + Jet_pt_4_ + len_*iCand) = tr.getVar<double>("Jet_pt_4");
        if(Jet_pt_5_ >= 0) *(basePtr_ + Jet_pt_5_ + len_*iCand) = tr.getVar<double>("Jet_pt_5");
        if(Jet_pt_6_ >= 0) *(basePtr_ + Jet_pt_6_ + len_*iCand) = tr.getVar<double>("Jet_pt_6");
        if(Jet_eta_1_ >= 0) *(basePtr_ + Jet_eta_1_ + len_*iCand) = tr.getVar<double>("Jet_eta_1");
        if(Jet_eta_2_ >= 0) *(basePtr_ + Jet_eta_2_ + len_*iCand) = tr.getVar<double>("Jet_eta_2");
        if(Jet_eta_3_ >= 0) *(basePtr_ + Jet_eta_3_ + len_*iCand) = tr.getVar<double>("Jet_eta_3");
        if(Jet_eta_4_ >= 0) *(basePtr_ + Jet_eta_4_ + len_*iCand) = tr.getVar<double>("Jet_eta_4");
        if(Jet_eta_5_ >= 0) *(basePtr_ + Jet_eta_5_ + len_*iCand) = tr.getVar<double>("Jet_eta_5");
        if(Jet_eta_6_ >= 0) *(basePtr_ + Jet_eta_6_ + len_*iCand) = tr.getVar<double>("Jet_eta_6");
        if(Jet_phi_1_ >= 0) *(basePtr_ + Jet_phi_1_ + len_*iCand) = tr.getVar<double>("Jet_phi_1");
        if(Jet_phi_2_ >= 0) *(basePtr_ + Jet_phi_2_ + len_*iCand) = tr.getVar<double>("Jet_phi_2");
        if(Jet_phi_3_ >= 0) *(basePtr_ + Jet_phi_3_ + len_*iCand) = tr.getVar<double>("Jet_phi_3");
        if(Jet_phi_4_ >= 0) *(basePtr_ + Jet_phi_4_ + len_*iCand) = tr.getVar<double>("Jet_phi_4");
        if(Jet_phi_5_ >= 0) *(basePtr_ + Jet_phi_5_ + len_*iCand) = tr.getVar<double>("Jet_phi_5");
        if(Jet_phi_6_ >= 0) *(basePtr_ + Jet_phi_6_ + len_*iCand) = tr.getVar<double>("Jet_phi_6");
        if(Jet_m_1_ >= 0) *(basePtr_ + Jet_m_1_ + len_*iCand) = tr.getVar<double>("Jet_m_1");
        if(Jet_m_2_ >= 0) *(basePtr_ + Jet_m_2_ + len_*iCand) = tr.getVar<double>("Jet_m_2");
        if(Jet_m_3_ >= 0) *(basePtr_ + Jet_m_3_ + len_*iCand) = tr.getVar<double>("Jet_m_3");
        if(Jet_m_4_ >= 0) *(basePtr_ + Jet_m_4_ + len_*iCand) = tr.getVar<double>("Jet_m_4");
        if(Jet_m_5_ >= 0) *(basePtr_ + Jet_m_5_ + len_*iCand) = tr.getVar<double>("Jet_m_5");
        if(Jet_m_6_ >= 0) *(basePtr_ + Jet_m_6_ + len_*iCand) = tr.getVar<double>("Jet_m_6");
        if(BestComboAvgMass_ >= 0) *(basePtr_ + BestComboAvgMass_ + len_*iCand) = tr.getVar<double>("BestComboAvgMass");
    }
};

class DeepEventShape
{
private:
    double discriminator_;
    std::string modelFile_, inputOp_, outputOp_;

    //Tensoflow session pointer
    TF_Session* session_;

    //Input variable names 
    std::vector<std::string> vars_;

    std::vector<TF_Output>     inputs_;
    std::vector<TF_Output>     outputs_;
    std::vector<TF_Operation*> targets_;

    //variable calclator
    std::shared_ptr<EventShapeCalculator> varCalculator_;

    void getParameters(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& localContextName)
    {
        //Construct contexts
        cfg::Context localCxt(localContextName);

        modelFile_     = cfgDoc->get("modelFile", localCxt, "");
        inputOp_       = cfgDoc->get("inputOp",   localCxt, "x");
        outputOp_      = cfgDoc->get("outputOp",  localCxt, "y");

        int iVar = 0;
        bool keepLooping;
        do
        {
            keepLooping = false;
        
            //Get variable name
            std::string varName = cfgDoc->get("mvaVar", iVar, localCxt, "");
        
            //if it is a non empty string save in vector
            if(varName.size() > 0)
            {
                keepLooping = true;
        
                vars_.push_back(varName);
                std::cout<<varName<<std::endl;
            }
            ++iVar;
        }
        while(keepLooping);
        
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
        varCalculator_.reset(new EventShapeCalculator());
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

        int iCand = 0;
        varCalculator_->calculateVars(tr, iCand);

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
        double discriminator = static_cast<double>(discriminators[iCand*TF_Dim(output_values[0], 1)]);

        for(auto tensor : input_values)  TF_DeleteTensor(tensor);
        for(auto tensor : output_values) TF_DeleteTensor(tensor);
        
        TF_DeleteStatus(status);

        //Define bins (try to have the same amount of events in all 4 bins)
        bool deepESM_bin1 = discriminator > 0.0000 && discriminator <= 0.4425;
        bool deepESM_bin2 = discriminator > 0.4425 && discriminator <= 0.4975;
        bool deepESM_bin3 = discriminator > 0.4975 && discriminator <= 0.5225;
        bool deepESM_bin4 = discriminator > 0.5225 && discriminator <= 1.0000;
        
        // Register Variables
        tr.registerDerivedVar("deepESM_val", discriminator);
        tr.registerDerivedVar("deepESM_bin1", deepESM_bin1);
        tr.registerDerivedVar("deepESM_bin2", deepESM_bin2);
        tr.registerDerivedVar("deepESM_bin3", deepESM_bin3);
        tr.registerDerivedVar("deepESM_bin4", deepESM_bin4);
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
    DeepEventShape(const std::string cfgFileName = "DeepEventShape.cfg", std::string localContextName = "Info")
    {
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

        std::cout<<"Using "+cfgFileName+" as the DeepEventShape config file"<<std::endl;
    }

    ~DeepEventShape()
    {
        //tensorflow status variable
        //TF_Status* status = TF_NewStatus();
        //TF_DeleteSession(session_, status);
        //TF_DeleteStatus(status);
    }
    
    void operator()(NTupleReader& tr)
    {
        runDeepEventShape(tr);
    }
};

#endif
