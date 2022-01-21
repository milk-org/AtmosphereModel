/**
 * @file    AtmosphereModel_Create_from_CONF.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_tools/COREMOD_tools.h"

#include "AtmosphereModel.h"
#include "AtmosphereModel_RefractionPath.h"
#include "AtmosphereModel_build_stdAtmModel.h"
#include "AtmosphereModel_stdAtmModel_ria.h"

//
// load refractive indices and abs coeff
// fname is the file name (eg. "RIA_O2.dat")
// lptr is pointer to wavelength
// RIptr is pointer to refractive index array
// absptr is pointer to absorption coeff array
//
long ATMOSPHEREMODEL_loadRIA(const char *restrict fname, RIAdata *riadat)
{
    FILE  *fp;
    long   nbpt = 0;
    double lmin, lmax;
    double v0, v1, v2;

    printf("Reading Refractive Index and Abs for %s\n", riadat->name);

    riadat->init = 0;
    riadat->NBpt = 0;
    riadat->Z    = 1.0; // default

    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("cannot open file \"%s\"\n", fname);
        return 0;
    }

    int r = fscanf(fp, "# %ld %lf %lf\n", &nbpt, &lmin, &lmax);
    if (r == 3)
    {
        printf("%ld points from %g m to %g m\n", nbpt, lmin, lmax);
    }
    else
    {
        printf("Read error\n");
        exit(0);
    }

    riadat->NBpt      = nbpt;
    riadat->lambdamin = lmin;
    riadat->lambdamax = lmax;

    if (riadat->NBpt > 0)
    {
        riadat->init   = 1;
        riadat->lambda = (double *) malloc(sizeof(double) * riadat->NBpt);
        riadat->rindex = (double *) malloc(sizeof(double) * riadat->NBpt);
        riadat->abs    = (double *) malloc(sizeof(double) * riadat->NBpt);

        for (long i = 0; i < riadat->NBpt; i++)
        {
            int r = fscanf(fp, "%lf %lf %lf\n", &v0, &v1, &v2);
            if (r == 3)
            {
                riadat->lambda[i] = v0;
                riadat->rindex[i] = v1;
                riadat->abs[i]    = v2;
            }
            else
            {
                printf("Read error\n");
                exit(0);
            }
        }
    }
    fclose(fp);

    return riadat->NBpt;
}

/// read configuration file and create atmosphere model

ATMOSPHERE_MODEL AtmosphereModel_Create_from_CONF(const char *restrict CONFFILE,
                                                  float slambda)
{
    ATMOSPHERE_MODEL atm;

    DEBUG_TRACEPOINT("Reading atmosphere parameters from file %s", CONFFILE);

    atm.ZenithAngle = read_config_parameter_float(CONFFILE, "ZENITH_ANGLE");
    atm.TimeDayOfYear =
        read_config_parameter_long(CONFFILE, "TIME_DAY_OF_YEAR");
    atm.TimeLocalSolarTime =
        read_config_parameter_long(CONFFILE, "TIME_LOCAL_SOLAR_TIME");

    atm.SiteLat  = read_config_parameter_float(CONFFILE, "SITE_LATITUDE");
    atm.SiteLong = read_config_parameter_float(CONFFILE, "SITE_LONGITUDE");
    atm.SiteAlt  = read_config_parameter_float(CONFFILE, "SITE_ALT");

    // temperature at site
    atm.SiteTPauto = read_config_parameter_float(CONFFILE, "SITE_TP_AUTO");

    // temperature at site
    atm.SiteTemp = read_config_parameter_float(CONFFILE, "SITE_TEMP");

    // pressure at site
    atm.SitePress = read_config_parameter_float(CONFFILE, "SITE_PRESS");

    atm.CO2_ppm = read_config_parameter_float(CONFFILE, "CO2_PPM");

    // water
    atm.SiteH2OMethod = read_config_parameter_long(CONFFILE, "SITE_H2O_METHOD");
    atm.SiteTPW       = read_config_parameter_float(CONFFILE, "SITE_TPW");
    atm.SiteRH        = read_config_parameter_float(CONFFILE, "SITE_RH");
    atm.SitePWSH      = read_config_parameter_float(CONFFILE, "SITE_PW_SCALEH");

    DEBUG_TRACEPOINT("Building standard atmoshpere model");
    AtmosphereModel_build_stdAtmModel(&atm, "atm.txt");

    printf("atm.SiteRH = %f\n", atm.SiteRH);

    atm.speciesRIA.NBspecies = atmNBspecies;

    DEBUG_TRACEPOINT("load refractive index data from %d species",
                     atm.speciesRIA.NBspecies);

    atm.speciesRIA.RIA_species =
        (RIAdata *) malloc(sizeof(RIAdata) * atm.speciesRIA.NBspecies);

    strcpy(atm.speciesRIA.RIA_species[speciesN2].name, "N2");
    strcpy(atm.speciesRIA.RIA_species[speciesO2].name, "O2");
    strcpy(atm.speciesRIA.RIA_species[speciesAr].name, "Ar");
    strcpy(atm.speciesRIA.RIA_species[speciesH2O].name, "H2O");
    strcpy(atm.speciesRIA.RIA_species[speciesCO2].name, "CO2");
    strcpy(atm.speciesRIA.RIA_species[speciesNe].name, "Ne");
    strcpy(atm.speciesRIA.RIA_species[speciesHe].name, "He");
    strcpy(atm.speciesRIA.RIA_species[speciesCH4].name, "CH4");
    strcpy(atm.speciesRIA.RIA_species[speciesKr].name, "Kr");
    strcpy(atm.speciesRIA.RIA_species[speciesH2].name, "H2");
    strcpy(atm.speciesRIA.RIA_species[speciesO3].name, "O3");
    strcpy(atm.speciesRIA.RIA_species[speciesN].name, "N");
    strcpy(atm.speciesRIA.RIA_species[speciesO].name, "O");
    strcpy(atm.speciesRIA.RIA_species[speciesH].name, "H");

    for (int sp = 0; sp < atm.speciesRIA.NBspecies; sp++)
    {
        char RIAfname[200];
        sprintf(RIAfname,
                "./RefractiveIndices/RIA_%s.dat",
                atm.speciesRIA.RIA_species[sp].name);
        ATMOSPHEREMODEL_loadRIA(RIAfname, &atm.speciesRIA.RIA_species[sp]);
    }

    DEBUG_TRACEPOINT("Set compressibility factors at STP");

    atm.speciesRIA.RIA_species[speciesN2].Z = 0.99971;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesCO2].Z = 0.99924;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesAr].Z = 0.99925;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesH2O].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesCO2].Z = 0.99435;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesNe].Z = 1.0005;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesHe].Z = 1.0005;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesCH4].Z = 0.99802;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesKr].Z = 0.99768;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesH2].Z = 1.0006;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesO3].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesN].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesO].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesH].Z = 1.0; // TO BE UPDATED

    if (0)
    {
        // ***************** refractive index as a function of lambda at site ****************************

        FILE *fp;

        DEBUG_TRACEPOINT("Write RindexSite.txt");
        fp = fopen("RindexSite.txt", "w");
        for (double lambda = 0.5e-6; lambda < 2e-6; lambda *= 1.0 + 1e-3)
        {
            RIAvalue ria =
                AtmosphereModel_stdAtmModel_ria(atm, atm.SiteAlt, lambda, 0);
            fprintf(fp, "%.8g %.14f %.14f\n", lambda, ria.rindex, ria.abscoeff);
            //printf("%10.6f    %.10f\n", lambda * 1.0e6, ria.rindex);
            fflush(stdout);
        }
        fclose(fp);
    }

    // *************** LIGHT PATH THROUGH ATMOSPHERE AT TWO SEPARATE WAVELENGTHS ***********************

    AtmosphereModel_RefractionPath(atm, 0.55e-6, atm.ZenithAngle, 1);
    EXECUTE_SYSTEM_COMMAND("mv refractpath.txt refractpath_0550.txt");

    AtmosphereModel_RefractionPath(atm, slambda, atm.ZenithAngle, 1);
    EXECUTE_SYSTEM_COMMAND("mv refractpath.txt refractpath_%04ld.txt",
                           (long) (1e9 * slambda + 0.5));

    if (0)
    {
        // refractive angle as a function of wavelength
        printf("Computing refraction angle ... \n");
        fflush(stdout);
        FILE *fp;
        fp = fopen("RefractAngle.dat", "w");
        for (double l = 0.5e-6; l < 10.0e-6; l *= 1.0 + 1e-2)
        {
            fprintf(fp,
                    "%20.18f  %.6f\n",
                    l,
                    AtmosphereModel_RefractionPath(atm, l, atm.ZenithAngle, 0));
        }
        fclose(fp);
        printf("done\n");
        fflush(stdout);
    }

    return atm;
}
