#ifndef DEEPEVENTSHAPE_H
#define DEEPEVENTSHAPE_H

#include "tensorflow/c/c_api.h"
#include "TopTagger/CfgParser/interface/Context.hh"
#include "TopTagger/CfgParser/interface/CfgDocument.hh"
#include "Framework/Framework/include/Utility.h"

#include "cstdlib"
#include "cstdio"
#include "cstring"

class EventShapeCalculator
{
private:
    float* basePtr_;
    std::vector<int> varIndex_;
    std::vector<std::string> vars_;
    std::string myVarSuffix_;
    std::vector<std::pair<std::string, int>> varPairs_;

public:
EventShapeCalculator(const std::vector<std::string>& vars, std::string myVarSuffix = "")
    : vars_(vars), myVarSuffix_(myVarSuffix)
    {
        for(const auto& v : vars_)
        {
            varPairs_.emplace_back(v, -1);
        }
    }
    /**
     *The job of mapVars is to populate the internal offests for all variables in the input variable list with their memory location in the data array.  To be called only once.
     */
    void mapVars(const std::vector<std::string>& vars)
    {
        for(unsigned int i = 0; i < vars.size(); ++i)
        {
            varPairs_[i].second = i;
        }
    }
    /**
     *The job of setPtr is to set the starting place of memory block where the data will be written. To be called only once for the creation of the array pointed to by data.
     */
    void setPtr(float* data) {basePtr_ = data;}
    /**
     *Calculate the requested variables and store the values directly in the input array for the MVA
     */
    template<typename T1, typename T2 = double> void calculateVar(const NTupleReader& tr, const int varId, const std::string& name)
    {
        if(varId >= 0) *(basePtr_ + varId) = static_cast<T1>( tr.getVar<T2>(name) );
    }

    void calculateVars(const NTupleReader& tr)
    {
        for(const auto& pair : varPairs_)
        {
            calculateVar<double>(tr, pair.second, pair.first+myVarSuffix_);
        }
    }
};

class DeepEventShape
{
private:
    double discriminator_;
    std::string modelFile_, inputOp_, outputOp_, year_, name_, nJetVar_, myVarSuffix_;
    int minNJet_, maxNJet_;
    bool firstEvent_;

    //Tensoflow session pointer
    TF_Session* session_;

    //Input variable names
    std::vector<std::string> outputOpVec_;
    std::vector<int> outputCmVec_;
    std::vector<std::string> vars_;
    std::map<std::string, std::vector<double> > binEdgesPerRegion_;
    std::vector<std::string> regions_;

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

        inputOp_     = cfgDoc->get("inputOp",   localCxt, "main_input");
        outputOp_    = cfgDoc->get("outputOp",  localCxt, "");
        outputOpVec_ = getVecFromCfg<std::string>(cfgDoc, "outputOpVec", localCxt, "");       
        outputCmVec_ = getVecFromCfg<int>(        cfgDoc, "outputCmVec", localCxt, -1);       
        year_        = cfgDoc->get("year",      localCxt, "");
        name_        = cfgDoc->get("name",      localCxt, "");
        nJetVar_     = cfgDoc->get("nJetVar",   localCxt, "NGoodJets");
        minNJet_     = cfgDoc->get("minNJet",   localCxt, 7);
        maxNJet_     = cfgDoc->get("maxNJet",   localCxt, 7);

        // Try looking for the array of input variables with name mvaVar[0], mvaVar[1], etc...
        vars_        = getVecFromCfg<std::string>(cfgDoc, "mvaVar", localCxt, "");
    
        // If no variables were found with array name mvaVar, then try inputVar
        if (vars_.size() == 0)
        {
            vars_ = getVecFromCfg<std::string>(cfgDoc, "inputVar", localCxt, "");
        }

        regions_     = getVecFromCfg<std::string>(cfgDoc, "regions", localCxt, "");

        for (const auto region : regions_) {
            binEdgesPerRegion_[region] = getVecFromCfg<double>(cfgDoc, "binEdges_"+region, localCxt, -1);
        }

        if(!outputOp_.empty())
        {
            outputOpVec_.push_back(outputOp_);
            outputCmVec_.push_back(2);
        }
        
        //Variable to hold tensorflow status
        TF_Status* status = TF_NewStatus();
        
        //get the grafdef from the file
        TF_Buffer* graph_def = read_file(modelFile_);
        
        // Import graph_def into graph
        TF_Graph* graph = TF_NewGraph();
        TF_ImportGraphDefOptions* graph_opts = TF_NewImportGraphDefOptions();
        TF_GraphImportGraphDef(graph, graph_def, graph_opts, status);
        TF_DeleteImportGraphDefOptions(graph_opts);
        if(TF_GetCode(status) != TF_OK) std::cerr<<utility::color("ERROR: Unable to import graph: "+std::string(TF_Message(status)), "red")<<std::endl;
        TF_DeleteBuffer(graph_def);
       
        //Create tensorflow session from imported graph
        TF_SessionOptions* sess_opts = TF_NewSessionOptions();
        uint8_t config[] = {0x10, 0x01};
        TF_SetConfig(sess_opts, static_cast<void*>(config), 2, status);
        session_ = TF_NewSession(graph, sess_opts, status);
        TF_DeleteSessionOptions(sess_opts);

        //Specify the name of the input layer
        TF_Operation* op_x  = TF_GraphOperationByName(graph, inputOp_.c_str());
        inputs_.emplace_back(TF_Output({op_x, 0}));

        //Specify the names of the output layers
        for(const auto& outName : outputOpVec_)
        {
            TF_Operation* op_y = TF_GraphOperationByName(graph, outName.c_str());
            outputs_.emplace_back(TF_Output({op_y, 0}));
            targets_.emplace_back(op_y);
        }

        //Clean up graph
        TF_DeleteGraph(graph);
        TF_DeleteStatus(status);

        //map variables
        varCalculator_.reset(new EventShapeCalculator(vars_, myVarSuffix_));
        varCalculator_->mapVars(vars_);
    }

    void runDeepEventShape(NTupleReader& tr)
    {
        //Check that the year the training is for is the same as file you are running over
        if(year_ != "" && firstEvent_)
        {
            const auto& runYear = tr.getVar<std::string>("runYear");
            try
            {                
                if(runYear != year_)
                {
                    throw "Warning: using DeepESM config file with \""+year_+"\" year but expected \""+runYear+"\" year";
                }
            }
            catch (const std::string msg) 
            {
                std::cerr<<utility::color(msg, "red")<<std::endl;
            }

            firstEvent_ = false;
        }

        auto& regions = tr.createDerivedVec<std::string>("regions_"+name_, regions_.size());
        for (unsigned int i = 0; i < regions_.size(); i++)
            regions.at(i) = regions_.at(i);

        //tensorflow status variable
        TF_Status* status = TF_NewStatus();
        
        //Create place to store the output vectors 
        std::vector<TF_Tensor*> output_values(outputs_.size());
       
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

        if(TF_GetCode(status) != TF_OK) std::cerr<<utility::color("ERROR: Unable to run graph: "+std::string(TF_Message(status)), "red")<<std::endl;

        //Get output discriminators 
        std::vector<std::vector<double>> discriminators(outputs_.size());

        for(unsigned int i = 0; i < output_values.size(); i++)
        {            
            auto* tensor = output_values[i];
            auto disc = static_cast<float*>(TF_TensorData(tensor));

            for(int j = 0; j < outputCmVec_[i]; j++)
            {
                discriminators[i].emplace_back(static_cast<double>(disc[j]));
            }

        }

        for(auto* tensor : input_values)  TF_DeleteTensor(tensor);
        for(auto* tensor : output_values) TF_DeleteTensor(tensor);        
        TF_DeleteStatus(status);

        // Register Variables and determine which regions/subregions an event is
        // based on inferencing with the Double DisCo NN
        double disc1   = discriminators[0][0];
        double disc2   = discriminators[0][2];
        double massReg = discriminators[1][0];
        tr.registerDerivedVar("DoubleDisCo_disc1_"+name_+myVarSuffix_, disc1);
        tr.registerDerivedVar("DoubleDisCo_disc2_"+name_+myVarSuffix_, disc2);
        tr.registerDerivedVar("DoubleDisCo_massReg_"+name_+myVarSuffix_, massReg);
        
        // Define and register deepESM bins
        const auto& NGoodJets = tr.getVar<int>(nJetVar_+myVarSuffix_);

        int kJet; // Iterator for edges, boundaries coming from Validation
        for (const auto region : regions_) {

            if (binEdgesPerRegion_.find(region) == binEdgesPerRegion_.end())
                continue;

            // Determine how many unique parameters per Njet bin per region
            int kSkip = binEdgesPerRegion_[region].size() / (maxNJet_-minNJet_+1);

            if      (NGoodJets <= minNJet_)
                kJet = 0;
            else if (NGoodJets <= maxNJet_)
                kJet = kSkip*(NGoodJets-minNJet_);
            else if (NGoodJets  > maxNJet_)
                kJet = kSkip*(maxNJet_-minNJet_);

            // Here, six parameters are extracted for each (Njets bin, ABCD) region in the order:
            // disc1, disc2, leftBoundary, rightBoundary, topBoundary, bottomBoundary
            double disc1edge      = binEdgesPerRegion_[region][kJet];
            double disc2edge      = binEdgesPerRegion_[region][kJet+1];
            double leftBoundary   = binEdgesPerRegion_[region][kJet+2];
            double rightBoundary  = binEdgesPerRegion_[region][kJet+3];
            double topBoundary    = binEdgesPerRegion_[region][kJet+4];
            double bottomBoundary = binEdgesPerRegion_[region][kJet+5];

            // Discriminant 1 placed on horizontal (left-to-right) axis
            // Discriminant 2 placed on vertical (top-to-bottom) axis

            // Check that the current event lives in the region in question.
            // Otherwise, it cannot possibly live in any of the "A", "B", "C", "D" subregions.
            bool withinBoundaries = disc1 >= leftBoundary   &&
                                    disc1 <= rightBoundary  && 
                                    disc2 >= bottomBoundary &&
                                    disc2 <= topBoundary;

            // Upper-right subregion
            bool passBinA = disc1 > disc1edge &&
                            disc2 > disc2edge &&
                            withinBoundaries;

            // Upper-left subregion
            bool passBinB = disc1 < disc1edge &&
                            disc2 > disc2edge &&
                            withinBoundaries;

            // Lower-right subregion
            bool passBinC = disc1 > disc1edge &&
                            disc2 < disc2edge &&
                            withinBoundaries;

            // Lower-left subregion
            bool passBinD = disc1 < disc1edge &&
                            disc2 < disc2edge &&
                            withinBoundaries;

            auto& passVec = tr.createDerivedVec<bool>("DoubleDisCo_"+region+"_"+name_+myVarSuffix_, 4);

            passVec.at(0) = passBinA;
            passVec.at(1) = passBinB;
            passVec.at(2) = passBinC;
            passVec.at(3) = passBinD;
        }
    }

    static void free_buffer(void* data, size_t) 
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
        , year_(husk.year_)
        , name_(husk.name_)
        , nJetVar_(husk.nJetVar_)
        , myVarSuffix_(husk.myVarSuffix_)
        , minNJet_(husk.minNJet_)
        , maxNJet_(husk.maxNJet_)
        , firstEvent_(husk.firstEvent_)
        , session_(husk.session_)
        , outputOpVec_(husk.outputOpVec_)
        , outputCmVec_(husk.outputCmVec_)
        , vars_(husk.vars_)
        , inputs_(husk.inputs_)
        , outputs_(husk.outputs_)
        , targets_(husk.targets_)
        , varCalculator_(husk.varCalculator_)
    {
        husk.session_ = nullptr;
    }
    
    DeepEventShape(const std::string& cfgFileName = "DeepEventShape.cfg", const std::string& modelFile = "keras_frozen.pb", const std::string& localContextName = "Info", 
                   const bool printStatus = true, const std::string& myVarSuffix = "")
        : modelFile_(modelFile)
        , myVarSuffix_(myVarSuffix)
        , firstEvent_(true)
    {
        if(printStatus) std::cout<<"Setting up DeepEventShape: Using "+cfgFileName+" and "+modelFile+" as the DeepEventShape config file and training file"<<std::endl;
        
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
