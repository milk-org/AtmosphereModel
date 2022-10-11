/**
 * @file    AirMixture_ria.c
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AtmosphereModel.h"
#include "OpticsMaterials/OpticsMaterials_n.h"

// densities in [cm-3]
//
// This routine assumes Z=1
//
RIAvalue
AirMixture_ria(ATM_SPECIES_RIADATA speciesRIA, double lambda, double *densarray)
{
    RIAvalue ria;

    // [m-3] for  1 atm (= 101.325 kPa) and 0 °C (= 273.15 K)
    double LoschmidtConstant = 2.6867805e25;

    double denstotal = 0.0;

    ria.abscoeff = 0.0;
    static int lliprecomp[atmNBspecies];

    double LLcumul = 0.0;

    for(int specindex = 0; specindex < atmNBspecies; specindex++)
    {
        double abscoeff = 0.0;
        double n;

        if(speciesRIA.RIA_species[specindex].init == 1)
        {
            long lli = lliprecomp[specindex];
            if(lli < 0)
            {
                lli = (long)(speciesRIA.RIA_species[specindex].NBpt / 2);
            }

            int llistep = 100;
            int llidir  = 1; // direction : -1=neg, 1=pos
            while(llistep != 1)
            {
                llistep = (long)(0.3 * llistep);
                if(llistep == 0)
                {
                    llistep = 1;
                }
                while(
                    (speciesRIA.RIA_species[specindex].lambda[lli] * llidir <
                     lambda * llidir) &&
                    (lli < speciesRIA.RIA_species[specindex].NBpt - llistep) &&
                    (lli > llistep))
                {
                    lli += llidir * llistep;
                }
                llidir = -llidir;
            }
            float alpha =
                (lambda - speciesRIA.RIA_species[specindex].lambda[lli]) /
                (speciesRIA.RIA_species[specindex].lambda[lli + 1] -
                 speciesRIA.RIA_species[specindex].lambda[lli]);
            //    printf("%d   alpha = %f \n", llidir, alpha);
            lliprecomp[specindex] = lli;
            n = (1.0 - alpha) * speciesRIA.RIA_species[specindex].rindex[lli] +
                alpha * speciesRIA.RIA_species[specindex].rindex[lli + 1];
            abscoeff = speciesRIA.RIA_species[specindex].abs[lli];
        }
        else
        {
            n = OpticsMaterials_n(
                    OpticsMaterials_code(speciesRIA.RIA_species[specindex].name),
                    lambda);
            //printf("    %20s  %.9f\n", speciesRIA.RIA_species[specindex].name, n);

            if(n < 0.0)  // not available
            {
                n = 1.0; // default
            }

            abscoeff = 0.0;
        }
        double LL = (n * n - 1) / (n * n + 2);
        double tmpc =
            densarray[specindex] /
            (LoschmidtConstant / 1e6 / speciesRIA.RIA_species[specindex].Z);
        LL *= tmpc;
        ria.abscoeff += tmpc * abscoeff;
        denstotal += densarray[specindex];
        LLcumul += LL;
        //printf("  %2d     %.6f  %g\n", specindex, LL, densarray[specindex]);
    }

    ria.rindex = sqrt((2.0 * LLcumul + 1.0) / (1.0 - LLcumul));

    //printf("ria.rindex = %.6f\n", ria.rindex);

    return ria;
}
