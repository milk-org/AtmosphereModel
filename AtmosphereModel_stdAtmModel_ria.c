/**
 * @file    AtmosphereModel_stdAtmModel_ria.c
 *
 */


#include <stdio.h>
#include <stdlib.h>

#include "CommandLineInterface/CLIcore.h"

#include "AtmosphereModel.h"
#include "AirMixture_ria.h"

//
// alt [m]
// lambda [m]
//
// computes n
// and absorption coefficient
//
//
RIAvalue AtmosphereModel_stdAtmModel_ria(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
)
{
    RIAvalue ria;
    
    float ifrac;

	// [m-3] for  1 atm (= 101.325 kPa) and 0 °C (= 273.15 K)
    double LoschmidtConstant = 2.6867805e25;


    
    double *densarray;

    long zindex = (long) (alt / ATM_VPROF_STEPSIZE);
    if(zindex > ATM_VPROF_NBSTEP-2)
    {
        zindex = ATM_VPROF_NBSTEP-1;
    }
    ifrac = 1.0 * alt / ATM_VPROF_STEPSIZE - zindex;

    if(ifrac < 0)
    {
        ifrac = 0.0;
    }

    if(ifrac > 1.0)
    {
        ifrac = 1.0;
    }
    
	DEBUG_TRACEPOINT("zindex = %ld", zindex);
    
    DEBUG_TRACEPOINT("atm.speciesRIA.NBspecies = %d", atm.speciesRIA.NBspecies);
    
    // densities [cm-3]
    double denstotal = 0.0;
    densarray = (double*) malloc(sizeof(double)*atm.speciesRIA.NBspecies);
    for(int spindex=0; spindex < atm.speciesRIA.NBspecies; spindex++)
    {
		densarray[spindex] = ((1.0 - ifrac) * atm.vprof_dens_species[spindex].val[zindex] + ifrac * atm.vprof_dens_species[spindex].val[zindex + 1]);
		denstotal += densarray[spindex];
	}
	
	DEBUG_TRACEPOINT("denstotal = %g", denstotal);

    ria = AirMixture_ria(atm.speciesRIA, lambda, densarray);

    if(mode == 1) // testing
    {
        printf("\n");
        printf("alt = %f m\n", alt);
        printf("zindex = %ld\n", zindex);
        printf("\n");
        printf("     N2  =  %5.3f\n", densarray[speciesN2] * 1.0e6 / LoschmidtConstant);
        printf("     O2  =  %5.3f\n", densarray[speciesO2] * 1.0e6 / LoschmidtConstant);
        printf("\n");
        printf("denstotal = %g    %g\n", denstotal,
               denstotal * 1.0e6 / LoschmidtConstant); // TEST
        printf("\n");
    }
    
    free(densarray);

    return ria;
}



double AtmosphereModel_stdAtmModel_N(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
)
{
	RIAvalue ria;
	
	ria = AtmosphereModel_stdAtmModel_ria(atm, alt, lambda, mode);
	
	return (1.0-ria.rindex);
}
