/**
 * @file    AtmosphereModel_save_stdAtmModel.c
 *
 */

#include <stdio.h>


#include "AtmosphereModel.h"


int AtmosphereModel_save_stdAtmModel(
    ATMOSPHERE_MODEL *atm,
    const char *restrict fname
)
{
    FILE *fp;

    fp = fopen(fname, "w");

    fprintf(fp,
            "#  1:alt[m]  2:denstot[part/cm3] 3:N2  4:O2  5:Ar  6:H2O  7:CO2  8:Ne  9:He  10:CH4  11:Kr  12:H2  13:N  14:O  15:H    16:density[g/cm3]  17:temperature[K] 18:pressure[stdatm]  19:RH\n");

    for(long zindex = 0; zindex < ATM_VPROF_NBSTEP; zindex++)
    {
        fprintf(fp,
                "%6.0f   %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g   %15.8g    %8.3lf   %12.10lf  %7.5f\n",
                ATM_VPROF_STEPSIZE * zindex,
                atm->vprof_denstot.val[zindex],
                atm->vprof_dens_species[0].val[zindex],
                atm->vprof_dens_species[1].val[zindex],
                atm->vprof_dens_species[2].val[zindex],
                atm->vprof_dens_species[3].val[zindex],
                atm->vprof_dens_species[4].val[zindex],
                atm->vprof_dens_species[5].val[zindex],
                atm->vprof_dens_species[6].val[zindex],
                atm->vprof_dens_species[7].val[zindex],
                atm->vprof_dens_species[8].val[zindex],
                atm->vprof_dens_species[9].val[zindex],
                atm->vprof_dens_species[10].val[zindex],
                atm->vprof_dens_species[11].val[zindex],
                atm->vprof_dens_species[12].val[zindex],
                atm->vprof_dens_species[13].val[zindex],
                atm->vprof_density.val[zindex],
                atm->vprof_temperature.val[zindex],
                atm->vprof_pressure.val[zindex],
                atm->vprof_RH.val[zindex]
               );
    }
    fclose(fp);

    return 0;
}
