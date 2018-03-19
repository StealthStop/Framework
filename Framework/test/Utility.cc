#include "Utility.h"
#include <math.h>

namespace utility
{
    double calcDPhi(double phi1, double phi2)
    {
        double dphi = phi1 - phi2 ;
        if ( dphi >  3.14159265 ) dphi -= 2*3.14159265 ;
        if ( dphi < -3.14159265 ) dphi += 2*3.14159265 ;
        return dphi;
    }
    
    double calcDR(double eta1, double eta2, double phi1, double phi2)
    {
        double deta = fabs( eta1 - eta2 ) ;
        
        double dphi = phi1 - phi2 ;
        if ( dphi > 3.1415926 ) dphi -= 2*3.1415926 ;
        if ( dphi <-3.1415926 ) dphi += 2*3.1415926 ;
        
        return sqrt( dphi*dphi + deta*deta ) ;
    }

    
}
