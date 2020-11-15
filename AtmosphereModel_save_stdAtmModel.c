/**
 * @file    AtmosphereModel_save_stdAtmModel.c
 *
 */

#include <stdio.h>


#include "AtmosphereModel.h"


int AtmosphereModel_save_stdAtmModel(
    ATMOSPHERE_MODEL atm,
    const char *restrict fname
)
{
    long i;
    FILE *fp;

    fp = fopen(fname, "w");

    fprintf(fp,
            "#  1:alt[m]  2:denstot[part/cm3] 3:N2  4:O2  5:Ar  6:H2O  7:CO2  8:Ne  9:He  10:CH4  11:Kr  12:H2  13:N  14:O  15:H    16:density[g/cm3]  17:temperature[K] 18:pressure[stdatm]  19:RH\n");

    for(i = 0; i < 10000; i++)
    {
        fprintf(fp,
                "%6.0f   %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g   %15.8g    %8.3lf   %12.10lf  %7.5f\n",
                10.0 * i,
                atm.vprof_denstot.val[i],
                atm.vprof_dens_species[0].val[i],
                atm.vprof_dens_species[1].val[i],
                atm.vprof_dens_species[2].val[i],
                atm.vprof_dens_species[3].val[i],
                atm.vprof_dens_species[4].val[i],
                atm.vprof_dens_species[5].val[i],
                atm.vprof_dens_species[6].val[i],
                atm.vprof_dens_species[7].val[i],
                atm.vprof_dens_species[8].val[i],
                atm.vprof_dens_species[9].val[i],
                atm.vprof_dens_species[10].val[i],
                atm.vprof_dens_species[11].val[i],
                atm.vprof_dens_species[12].val[i],
                atm.vprof_dens_species[13].val[i],
                atm.vprof_density.val[i],
                atm.vprof_temperature.val[i],
                atm.vprof_pressure.val[i],
                atm.vprof_RH.val[i]
               );
    }
    fclose(fp);

    return 0;
}
