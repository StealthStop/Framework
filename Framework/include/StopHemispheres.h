#ifndef StopHemispheres_h
#define StopHemispheres_h

#include <vector>
#include <iostream>
#include <cmath>

class StopHemispheres 
{

// There are 2 constructors:
// 1. Constructor 
//    taking as argument vectors of Px, Py, Pz and E of the objects in the event that should be separated, 
//    the seeding method and the hemisphere association method,
// 2. Constructor 
//    taking as argument vectors of Px, Py, Pz and E of the objects in the event that should be separated. 
//    the seeding method and the hemisphere association method should then be defined by SetMethod(seeding_method, association_method).
// 
// Seeding method to choice of 2 inital axes:
//   2: 2 objects who give maximal invariant mass (recommended)
//
// StopHemispheres association method:
//   3: minimal Lund distance (recommended)
//
// Note that SetMethod() also allows the seeding and/or association method to be redefined for an existing hemisphere object.
// The GetAxis() or GetGrouping() is  then recomputed using the newly defined methods.
//
// Some parameters possibly affecting the logic of hemisphere reconstruction can also be set by the user:
//   SetNoSeed()      : prevent the given object from being used as a seed (but it can be associated to hemispheres).
//   SetNoAssoc()     : prevent the given object from being used in the hemisphere association (then it cannot be used as seed either).
//   ClearAllNoLists(): reset the list of NoSeed and NoAssoc objects to empty.
//   SetnItermax()    : maximum number of iterations allowed for the association methods (default = 100).
    
public:
    enum SeedMethod {NoSeed = 0, InvMassSeed = 2, TopSeed = 5};

    StopHemispheres(std::vector<float> Px_vector, std::vector<float> Py_vector, std::vector<float> Pz_vector, std::vector<float> E_vector, SeedMethod seed_method, int hemisphere_association_method);
    StopHemispheres(std::vector<float> Px_vector, std::vector<float> Py_vector, std::vector<float> Pz_vector, std::vector<float> E_vector);

    ~StopHemispheres(){};

    // where Nx, Ny, Nz are the direction cosines e.g. Nx = Px/P; P = momentum; E = energy
    // return Nx, Ny, Nz, P, E of the axis of group 1
    std::vector<float> getAxis1(); 
    // return Nx, Ny, Nz, P, E of the axis of group 2
    std::vector<float> getAxis2(); 

    // return vector with "1" and "2"'s according to which group the object belongs 
    // and 0 if the object is not associated to any hemisphere
    // (order of objects in vector is same as input)
    std::vector<int> getGrouping();

    // set or overwrite the seed and association methods
    void SetMethod(SeedMethod seed_method, int hemisphere_association_method) 
    {
        seed_meth = seed_method;
        hemi_meth = hemisphere_association_method;
        status    = 0;
    }

    // prevent an object from being used as a seed (but it can be associated to hemispheres)
    void SetNoSeed(int object_number) 
    {
        Object_Noseed[object_number] = 1; 
        status                       = 0;
    } 

    // prevent an object from being used for hemisphere asociation (and for seeding)
    void SetNoAssoc(int object_number) 
    {
        Object_Noassoc[object_number] = 1; 
        Object_Noseed[object_number]  = 1; 
        status                        = 0;
    } 

    // reset the list of NoSeed and NoAssoc objects to empty
    void ClearAllNoLists() 
    {
        for (int i = 0; i < (int)Object_Noseed.size(); ++i) 
        {
            Object_Noassoc[i] = 0; 
            Object_Noseed[i]  = 0; 
            status            = 0;
        }
    }

    // set the maximum number of iterations allowed for association
    void SetnItermax(int niter) 
    {
        nItermax = niter;
    }

    // controls the level of debug prints
    void SetDebug(int debug) { dbg = debug; } 
    int  GetNumLoop()        {return numLoop; }

private:
    int Reconstruct();

    std::vector<float> Object_Px;
    std::vector<float> Object_Py;
    std::vector<float> Object_Pz;
    std::vector<float> Object_P;
    std::vector<float> Object_Pt;
    std::vector<float> Object_E;
    std::vector<float> Object_Phi;
    std::vector<float> Object_Eta;
    std::vector<int> Object_Group;
    std::vector<int> Object_Noseed;
    std::vector<int> Object_Noassoc;

    std::vector<float> Axis1;
    std::vector<float> Axis2;

    SeedMethod seed_meth;
    int hemi_meth;
    int status;
    int nItermax;
    int dbg;
    int numLoop;
};

#endif
