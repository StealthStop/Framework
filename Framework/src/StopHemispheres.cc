#include "Framework/Framework/include/StopHemispheres.h"
#include "Framework/Framework/include/HemispheresUtilities.h"
#include <math.h>

using namespace std;

using std::vector;
using std::cout;
using std::endl;

// constructor specifying the seed and association methods
StopHemispheres::StopHemispheres(vector<float> Px_vector, vector<float> Py_vector, vector<float> Pz_vector,
                       vector<float> E_vector, SeedMethod seed_method, int hemisphere_association_method) : 
                       Object_Px(Px_vector), Object_Py(Py_vector), Object_Pz(Pz_vector), Object_E(E_vector), 
                       seed_meth(seed_method), hemi_meth(hemisphere_association_method), 
                       status(0), nItermax(100), dbg(0)  
{	
    for(int i = 0; i < (int) Object_Px.size(); i++)
    {
        Object_Noseed.push_back(0);
        Object_Noassoc.push_back(0);
    }
    numLoop =0;
}

// constructor without specification of the seed and association methods
// in this case, the latter must be given by calling SetMethod before invoking reconstruct()
StopHemispheres::StopHemispheres(vector<float> Px_vector, vector<float> Py_vector, vector<float> Pz_vector, vector<float> E_vector) : 
                       Object_Px(Px_vector), Object_Py(Py_vector), Object_Pz(Pz_vector), Object_E(E_vector), 
                       seed_meth(NoSeed), hemi_meth(0), 
                       status(0), nItermax(100), dbg(0)  
{	
    for(int i = 0; i < (int) Object_Px.size(); i++)
    {
		Object_Noseed.push_back(0);
		Object_Noassoc.push_back(0);
    }
    numLoop =0;
}

vector<float> StopHemispheres::getAxis1()
{
    if (status != 1) 
    {
        this->Reconstruct();
    }
    return Axis1;
}

vector<float> StopHemispheres::getAxis2()
{
    if (status != 1) 
    {
        this->Reconstruct();
    }
    return Axis2;
}

vector<int> StopHemispheres::getGrouping()
{
    if (status != 1) 
    {
        this->Reconstruct();
    }
    return Object_Group;
}

int StopHemispheres::Reconstruct()
{

    // definition of the vectors used internally:
    // Object_xxx                 : xxx = Px, Py, Pz, E for input values
    //                              xxx = P, Pt, Eta, Phi, Group for internal use
    // Axis1                      : final hemisphere axis 1
    // Axis2                      : final hemisphere axis 2
    // Sum1_xxx                   : hemisphere 1 being updated during the association iterations
    // Sum2_xxx                   : hemisphere 2 being updated during the association iterations
    // NewAxis1_xxx, NewAxis2_xxx : temporary axes for calculation in association methods 2 and 3

    // initialize numLoop for Zero
    numLoop = 0;

    int vsize = (int) Object_Px.size();
    if((int) Object_Py.size() != vsize || (int) Object_Pz.size() != vsize)
    {
        cout << "WARNING!!!!! Input vectors have different size! Fix it!" << endl;
        return 0;
    }

    if (dbg > 0) 
    {
        cout << " StopHemispheres method " << endl;
    }

    // clear some vectors if method reconstruct() is called again
    if (!Object_P.empty())
    {
        Object_P.clear();
        Object_Pt.clear();
        Object_Eta.clear();
        Object_Phi.clear();
        Object_Group.clear();
        Axis1.clear();
        Axis2.clear();
    }

    // initialize the vectors
    for(int j = 0; j < vsize; ++j)
    {
        Object_P.push_back(0);
        Object_Pt.push_back(0);
        Object_Eta.push_back(0);
        Object_Phi.push_back(0);
        Object_Group.push_back(0);
    }
	
    for(int j = 0; j < 5; ++j)
    {
        Axis1.push_back(0);
        Axis2.push_back(0);
    }
    
    // compute additional quantities for vectors Object_xxx
    float theta;
    for (int i = 0; i <vsize; ++i)
    {
        Object_P[i] = sqrt( (Object_Px[i] * Object_Px[i]) + (Object_Py[i] * Object_Py[i]) + (Object_Pz[i] * Object_Pz[i]) );
        if (Object_P[i] > Object_E[i]+0.001) 
        {
            cout << "WARNING!!!!! Object " << i << " has E = " << Object_E[i] << " less than P = " << Object_P[i] << " *** Fix it!" << endl;
            return 0;
        }
        
        Object_Pt[i] = sqrt( (Object_Px[i] * Object_Px[i]) + (Object_Py[i] * Object_Py[i]) );
        if (fabs(Object_Pz[i]) > 0.001) 
        {
            theta = atan(sqrt( (Object_Px[i] * Object_Px[i]) + (Object_Py[i] * Object_Py[i]) ) / Object_Pz[i]);
        } 
        else 
        {
            theta = M_PI/2;
        }
		
        if (theta < 0.) 
        { 
            theta = theta + M_PI; 
        }

        Object_Eta[i] = -log(tan(0.5*theta));
        Object_Phi[i] = atan2(Object_Py[i], Object_Px[i]);
        if (dbg > 0) 
        {
            cout << " Object " << i << " Eta = " << Object_Eta[i] << " Phi = " << Object_Phi[i] << endl;
        }
    }

    if (dbg > 0) 
    {
        cout << endl;
        cout << " Seeding method = " << seed_meth << endl;
    }
	
    // ------------------------------------------------------------------------------------
    // (1) SEED METHODS 
    // ------------------------------------------------------------------------------------
    // I_Max and J_Max are indices of the seeds in the vectors
    int I_Max = -1;
    int J_Max = -1;

    // ------------------------------------------- 
	// -- determine the seeds for seed method 2 
	//    -- OldSeed
	// ------------------------------------------- 
    if (seed_meth == InvMassSeed) 
    {
        float Mass_Max = 0.;
        float InvariantMass = 0.;
        
        // maximize the invariant mass of two objects
        for (int i = 0; i < vsize; ++i)
        {
            Object_Group[i] = 0;
            if (Object_Noseed[i] == 0)  
            {
                for (int j = i+1; j < vsize; ++j)
                {
                    if (Object_Noseed[j] == 0)
                    {
                        // either the invariant mass
                        if (seed_meth == InvMassSeed)
                        {
                            InvariantMass = (Object_E[i] + Object_E[j]) * (Object_E[i] + Object_E[j])
                                - (Object_Px[i] + Object_Px[j]) * (Object_Px[i] + Object_Px[j])
                                - (Object_Py[i] + Object_Py[j]) * (Object_Py[i] + Object_Py[j])
                                - (Object_Pz[i] + Object_Pz[j]) * (Object_Pz[i] + Object_Pz[j]);
                        }
                        
                        if (Mass_Max < InvariantMass)
                        {
                            Mass_Max = InvariantMass;
                            I_Max = i;
                            J_Max = j;
                        }
                    }
                }
            }
        }

        // if both seeds are found, save them as initial hemisphere axes
        if (J_Max > 0) 
        {
            Axis1[0] = Object_Px[I_Max] /  Object_P[I_Max];
            Axis1[1] = Object_Py[I_Max] /  Object_P[I_Max];
            Axis1[2] = Object_Pz[I_Max] /  Object_P[I_Max];
            Axis1[3] = Object_P[I_Max];
            Axis1[4] = Object_E[I_Max];

            Axis2[0] = Object_Px[J_Max] /  Object_P[J_Max];
            Axis2[1] = Object_Py[J_Max] /  Object_P[J_Max];
            Axis2[2] = Object_Pz[J_Max] /  Object_P[J_Max];
            Axis2[3] = Object_P[J_Max];
            Axis2[4] = Object_E[J_Max];
        } else 
        {
            return 0;
        }
		
        if (dbg > 0) 
        {
            cout << " Axis 1 is Object = " << I_Max << endl;
            cout << " Axis 2 is Object = " << J_Max << endl;
        }
    } 

    // -------------------------------------------
    // -- determine the seeds for seed method 5
    //    -- TopSeed
    // -------------------------------------------  
    else if ( seed_meth == TopSeed ) 
    {
        I_Max = 0;
        J_Max = 1;
        
        if (Object_Pt.size() > 1)
        {
            if ( Object_Pt[J_Max] > Object_Pt[I_Max] )
            {
                std::swap(I_Max, J_Max);
            }

            Axis1[0] = Object_Px[I_Max] / Object_P[I_Max];
            Axis1[1] = Object_Py[I_Max] / Object_P[I_Max];
            Axis1[2] = Object_Pz[I_Max] / Object_P[I_Max];
       	    Axis1[3] = Object_P[I_Max];
       	    Axis1[4] = Object_E[I_Max];

            Axis2[0] = Object_Px[J_Max] / Object_P[J_Max];
            Axis2[1] = Object_Py[J_Max] / Object_P[J_Max];
            Axis2[2] = Object_Pz[J_Max] / Object_P[J_Max];
       	    Axis2[3] = Object_P[J_Max];
       	    Axis2[4] = Object_E[J_Max];
        } else
        {
            cout << "Can not get 2 objects for TopSeed method !!!" << endl;
            return 0;
        }

        if ( dbg > 0 )
        {
            cout << " Axis 1 is Object = " << I_Max  << " with Pt " << Object_Pt[I_Max]<< endl;
            cout << " Axis 2 is Object = " << J_Max  << " with Pt " << Object_Pt[J_Max]<< endl;
        }

    } else if ( !(seed_meth == NoSeed) ) 
    {
        cout << "Please give a valid seeding method!" << endl;
        return 0;
    }
	    
    // ------------------------------------------------------------------------------------
    // (2) HEMISPHERE ASSOCIATION METHOD
    // ------------------------------------------------------------------------------------
    if (dbg > 0) 
    {
        cout << endl;
        cout << " Association method = " << hemi_meth << endl;
    }
    
    // ---------------------------------------------------------------------------
    // -- iterate to associate all objects to hemispheres (methods 1 to 3 only) 
    //    -- until no objects are moved from one to the other hemisphere 
    //    -- or the maximum number of iterations is reached
    // ---------------------------------------------------------------------------
    bool I_Move = true;

    while (I_Move && (numLoop < nItermax))
    {
        I_Move = false;
        numLoop++;
        if (dbg > 0) 
        {
            cout << " Iteration = " << numLoop << endl;
        }
        if(numLoop == nItermax-1)
        {
            cout << " Hemishpere: warning - reaching max number of iterations " << endl;
        }

        // initialize the current sums of Px, Py, Pz, E for the two hemispheres
        float Sum1_Px = 0.;
        float Sum1_Py = 0.;
        float Sum1_Pz = 0.;
        float Sum1_E  = 0.;
        float Sum2_Px = 0.;
        float Sum2_Py = 0.;
        float Sum2_Pz = 0.;
        float Sum2_E  = 0.;

        // -----------------------------------------------
	    // -- associate the objects for Lund distance
	    //    -- (method 3) 
	    // -----------------------------------------------
        if (hemi_meth == 3) 
        {            
            for (int i = 0; i < vsize; ++i)
            {
                // add the seeds to the sums, as they remain fixed
                if (i == I_Max) 
                {
                    Object_Group[i] = 1;
                    Sum1_Px         += Object_Px[i];
                    Sum1_Py         += Object_Py[i];
                    Sum1_Pz         += Object_Pz[i];
                    Sum1_E          += Object_E[i];
                } 
                
                else if (i == J_Max) 
                {
                    Object_Group[i] = 2;
                    Sum2_Px         += Object_Px[i];
                    Sum2_Py         += Object_Py[i];
                    Sum2_Pz         += Object_Pz[i];
                    Sum2_E          += Object_E[i];
                }
                // ------
                // sum 
                // ------
                else 
                {
                    if (Object_Noassoc[i] == 0)
                    {
                        // only 1 object maximum is moved in a given iteration
                        if(!I_Move)
                        {
                            // initialize the new hemispheres as the current ones
                            float NewAxis1_Px = Axis1[0] * Axis1[3];
                            float NewAxis1_Py = Axis1[1] * Axis1[3];
                            float NewAxis1_Pz = Axis1[2] * Axis1[3];
                            float NewAxis1_E  = Axis1[4];
                            float NewAxis2_Px = Axis2[0] * Axis2[3];
                            float NewAxis2_Py = Axis2[1] * Axis2[3];
                            float NewAxis2_Pz = Axis2[2] * Axis2[3];
                            float NewAxis2_E  = Axis2[4];
				            
                            // subtract the object from its hemisphere
                            if (Object_Group[i] == 1)
                            {
                                NewAxis1_Px = NewAxis1_Px - Object_Px[i];
                                NewAxis1_Py = NewAxis1_Py - Object_Py[i];
                                NewAxis1_Pz = NewAxis1_Pz - Object_Pz[i];
                                NewAxis1_E  = NewAxis1_E  - Object_E[i];
                            } 
                            else if (Object_Group[i] == 2) 
                            {
                                NewAxis2_Px = NewAxis2_Px - Object_Px[i];
                                NewAxis2_Py = NewAxis2_Py - Object_Py[i];
                                NewAxis2_Pz = NewAxis2_Pz - Object_Pz[i];
                                NewAxis2_E  = NewAxis2_E  - Object_E[i];
                            }
                           
                            // compute the invariant mass squared with each hemisphere
                            float mass1 =  NewAxis1_E - ( ( (Object_Px[i] * NewAxis1_Px) + (Object_Py[i] * NewAxis1_Py) + (Object_Pz[i] * NewAxis1_Pz) ) / Object_P[i] );
                            float mass2 =  NewAxis2_E - ( ( (Object_Px[i] * NewAxis2_Px) + (Object_Py[i] * NewAxis2_Py) + (Object_Pz[i] * NewAxis2_Pz) ) / Object_P[i] );

                            // ------------------------------
                            // -- Lund distance (method 3)
                            // ------------------------------
                            if (hemi_meth == 3) 
                            {
                                mass1 *= NewAxis1_E / ( (NewAxis1_E + Object_E[i]) * (NewAxis1_E + Object_E[i]) );
                                mass2 *= NewAxis2_E / ( (NewAxis2_E + Object_E[i]) * (NewAxis1_E + Object_E[i]) );
                            }

                            // and associate the object to the best hemisphere and add it to the sum
                            if (mass1 < mass2) 
                            {
                                if (Object_Group[i] != 1)
                                {
                                    I_Move = true;
                                }
                                Object_Group[i] = 1;
                                Sum1_Px         += Object_Px[i];
                                Sum1_Py         += Object_Py[i];
                                Sum1_Pz         += Object_Pz[i];
                                Sum1_E          += Object_E[i];
                            } 
                            else 
                            {
                                if (Object_Group[i] != 2)
                                {
                                    I_Move = true;
                                }
                                Object_Group[i] = 2;
                                Sum2_Px         += Object_Px[i];
                                Sum2_Py         += Object_Py[i];
                                Sum2_Pz         += Object_Pz[i];
                                Sum2_E          += Object_E[i];
                            }
                        } else 
                        {
                            if (Object_Group[i] == 1)
                            {
                                Sum1_Px += Object_Px[i];
                                Sum1_Py += Object_Py[i];
                                Sum1_Pz += Object_Pz[i];
                                Sum1_E  += Object_E[i];
                            } 
                            else if (Object_Group[i] == 2)
                            {
                                Sum2_Px += Object_Px[i];
                                Sum2_Py += Object_Py[i];
                                Sum2_Pz += Object_Pz[i];
                                Sum2_E  += Object_E[i];
                            }
                        }
                    }
                } // end loop over objects, Sum1_ and Sum2_ are now the updated hemispheres
            } // for loop
        } // hemisphere method 3
 
        else 
        {
            cout << "Please give a valid hemisphere association method!" << endl;
            return 0;
        }

        // ---------------------------------------------    	
        // -- recomputing the axes for next iteration
        // --------------------------------------------- 
        Axis1[3] = sqrt( (Sum1_Px * Sum1_Px) + (Sum1_Py * Sum1_Py) + (Sum1_Pz * Sum1_Pz) );
        if (Axis1[3] < 0.0001) 
        {
            cout << "ZERO objects in group 1! " << endl;
        } else 
        {
            Axis1[0] = Sum1_Px / Axis1[3];
            Axis1[1] = Sum1_Py / Axis1[3];
            Axis1[2] = Sum1_Pz / Axis1[3];
            Axis1[4] = Sum1_E;
        }

        Axis2[3] = sqrt( (Sum2_Px * Sum2_Px) + (Sum2_Py * Sum2_Py) + (Sum2_Pz * Sum2_Pz) );
        if (Axis2[3] < 0.0001) 
        {
            cout << " ZERO objects in group 2! " << endl;
        } else 
        {
            Axis2[0] = Sum2_Px / Axis2[3];
            Axis2[1] = Sum2_Py / Axis2[3];
            Axis2[2] = Sum2_Pz / Axis2[3];
            Axis2[4] = Sum2_E;
        }

        if (dbg > 0) 
        {
            cout << " Grouping = ";
            for (int i=0; i < vsize; i++)
            {
                cout << "  " << Object_Group[i];
            }
            cout << endl;
        }

    } // end of iteration
    
    status = 1;
    return 1;
}

