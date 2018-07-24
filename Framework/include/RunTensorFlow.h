#ifndef RUNTENSORFLOW_H
#define RUNTENSORFLOW_H

#include "tensorflow/c/c_api.h"
#include "TopTagger/CfgParser/include/Context.hh"
#include "TopTagger/CfgParser/include/CfgDocument.hh"

class RunTensorFlow
{
private:
    double discriminator_;
    double discOffset_;
    double discSlope_;
    std::string modelFile_, inputOp_, outputOp_;
    double csvThreshold_;
    double bEtaCut_;
    int maxNbInTop_;
    int NConstituents_;

    //Tensoflow session pointer
    TF_Session* session_;

    //Input variable names 
    std::vector<std::string> vars_;

    std::vector<TF_Output>     inputs_;
    std::vector<TF_Output>     outputs_;
    std::vector<TF_Operation*> targets_;

    void getParameters(const std::unique_ptr<cfg::CfgDocument>& cfgDoc, const std::string& localContextName)
    {
        //Construct contexts
        cfg::Context commonCxt("Common");
        cfg::Context localCxt(localContextName);

        discriminator_ = cfgDoc->get("discCut",       localCxt, -999.9);
        discOffset_    = cfgDoc->get("discOffset",    localCxt, 999.9);
        discSlope_     = cfgDoc->get("discSlope",     localCxt, 0.0);
        modelFile_     = cfgDoc->get("modelFile",     localCxt, "");
        inputOp_       = cfgDoc->get("inputOp",       localCxt, "x");
        outputOp_      = cfgDoc->get("outputOp",      localCxt, "y");
        NConstituents_ = cfgDoc->get("NConstituents", localCxt, 3);
        
        csvThreshold_  = cfgDoc->get("csvThreshold", localCxt, -999.9);
        bEtaCut_       = cfgDoc->get("bEtaCut",      localCxt, -999.9);
        maxNbInTop_    = cfgDoc->get("maxNbInTop",   localCxt, -1);

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
    }

    void runTensorFlow(NTupleReader& tr)
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
        //varCalculator_->setPtr(static_cast<float*>(TF_TensorData(input_values_0)));

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
        int iCand = 0;
        
        //discriminators is a 2D array, we only want the first entry of every array
        double discriminator = static_cast<double>(discriminators[iCand*TF_Dim(output_values[0], 1)]);

        for(auto tensor : input_values)  TF_DeleteTensor(tensor);
        for(auto tensor : output_values) TF_DeleteTensor(tensor);

        TF_DeleteStatus(status);

        // Register Variables
        tr.registerDerivedVar("tensorflow_val", discriminator);
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
    RunTensorFlow()
    {
        std::string cfgText = "";
        std::string localContextName = "";
        std::unique_ptr<cfg::CfgDocument> cfgDoc = cfg::CfgDocument::parseDocument(cfgText);
        getParameters(cfgDoc, localContextName);
    }

    ~RunTensorFlow()
    {
        //tensorflow status variable
        TF_Status* status = TF_NewStatus();
        TF_DeleteSession(session_, status);
        TF_DeleteStatus(status);
    }
    
    void operator()(NTupleReader& tr)
    {
        runTensorFlow(tr);
    }
};

#endif
