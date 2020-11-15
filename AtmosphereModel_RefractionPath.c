/**
 * @file    AtmosphereModel_RefractionPath.c
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "AtmosphereModel.h"
#include "AtmosphereModel_stdAtmModel_ria.h"


/// lambda in m
/// if WritePath = 1, write light path in output ASCII file
/// return atm refraction value [arcsec]
double AtmosphereModel_RefractionPath(
	ATMOSPHERE_MODEL atm,
    double lambda,
    double Zangle,
    int WritePath
)
{
    double h0, h1; // elevation (from Earth center direction, 0 at Earth surface)
    double lstep = 20.0; // beam distance step [m]
    double Re = 6371000.0; // Earth radius [m]
    double alphae; // angle to center of Earth
    double x0, y0, x1, y1; // cartesian coordinates
    double n0, n1;
    double alpha0, alpha1;
    double alpha; // current zenith angle
  
    FILE *fp;
    double pathl;
    double offsetangle;
    double Zangle0;
    long iter = 0;
    double errV = 100000.0;
    double flux;


    Zangle0 = Zangle;
    offsetangle = 0.0;

    while((iter < 20) && (errV > 0.00001))
    {
        //printf("----------------- iter %ld ---------------\n", iter);
        // upward path
        x0 = 0.0;
        y0 = 0.0;
        alpha = Zangle0;
        h0 = atm.SiteAlt;
        alpha1 = 0.0;
        RIAvalue riav0 = AtmosphereModel_stdAtmModel_ria(atm, h0, lambda, 0);
        n0 = riav0.rindex;

        if(WritePath)
        {
            fp = fopen("refractpath.txt", "w");
        }
        pathl = 0.0;
        h1 = 0.0;
        flux = 1.0;
        while(h1 < 99000.0)
        {
            x1 = x0 + lstep * sin(alpha);
            y1 = y0 + lstep * cos(alpha);

            h1 = sqrt(x1 * x1 + (y1 + atm.SiteAlt + Re) * (y1 + atm.SiteAlt + Re)) - Re;
            RIAvalue riav = AtmosphereModel_stdAtmModel_ria(atm, h1, lambda, 0);
            n1 = riav.rindex;

            alphae = atan2(x1, y1 + atm.SiteAlt + Re);
            flux *= exp(-lstep * riav.abscoeff);
            //printf("[%f %f] abscoeff =  %g      flux -> %lf\n", h1, lambda, v_ABSCOEFF, flux);
            //exit(0);
            alpha0 = alpha - alphae;
            alpha1 = asin(n0 * sin(alpha0) / n1);
            alpha = alpha1 + alphae;

            n0 = n1;
            x0 = x1;
            y0 = y1;
            h0 = h1;

            pathl += lstep;
            if(WritePath)
            {
                fprintf(fp, "%10.3f %10.3f %10.3f  %12.9f  %20g    %20g   %20g\n", h0, x0, y0,
                        alpha, alpha - Zangle, (alpha - Zangle) / M_PI * 180.0 * 3600.0,
                        y0 * tan(Zangle) - x0);
            }
        }
        if(WritePath)
        {
            fclose(fp);
        }
        offsetangle = alpha - Zangle;
        //printf("refraction error = %g arcsec\n", offsetangle/M_PI*180.0*3600.0);
        //printf("path = %f\n", pathl);
        Zangle0 -= offsetangle;
        errV = fabs(offsetangle / M_PI * 180.0 * 3600.0);
        iter++;
    }
    
    
    RIAvalue riav = AtmosphereModel_stdAtmModel_ria(atm, atm.SiteAlt, lambda, 0);
    
    
    printf("%10.6f um   %12.10f     Atmospheric Refraction = %10.6f arcsec   ",
           lambda * 1e6, riav.rindex,
           (Zangle - Zangle0) / M_PI * 180.0 * 3600.0);
    printf("TRANSMISSION = %lf\n", flux);
    //v_TRANSM = flux;

    return((Zangle - Zangle0) / M_PI * 180.0 * 3600.0);
}

