/**
 * @file    AtmosphereModel_Create_from_CONF.c
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "AtmosphereModel.h"
#include "AtmosphereModel_build_stdAtmModel.h"
#include "AtmosphereModel_stdAtmModel_ria.h"
#include "AtmosphereModel_RefractionPath.h"





//
// load refractive indices and abs coeff
// fname is the file name (eg. "RIA_O2.dat")
// lptr is pointer to wavelength
// RIptr is pointer to refractive index array
// absptr is pointer to absorption coeff array
//
long ATMOSPHEREMODEL_loadRIA(
    const char *restrict fname,
    RIAdata *riadat
)
{
    FILE *fp;
    long nbpt = 0;
    double lmin, lmax;
    double v0, v1, v2;

    printf("Reading Refractive Index and Abs for %s\n", riadat->name);

	riadat->init = 0;
	riadat->NBpt = 0;
	riadat->Z = 1.0; // default
	
    fp = fopen(fname, "r");
    if(fp == NULL)
    {
        printf("cannot open file \"%s\"\n", fname);
        return 0;
    }

    int r = fscanf(fp, "# %ld %lf %lf\n", &nbpt, &lmin, &lmax);
    if(r == 3)
    {
        printf("%ld points from %g m to %g m\n", nbpt, lmin, lmax);
    }
    else
    {
        printf("Read error\n");
        exit(0);
    }

	riadat->NBpt = nbpt;
	riadat->lambdamin = lmin;
	riadat->lambdamax = lmax;
	
    if(riadat->NBpt > 0)
    {
		riadat->init = 1;
        riadat->lambda = (double *) malloc(sizeof(double) * riadat->NBpt);
        riadat->rindex = (double *) malloc(sizeof(double) * riadat->NBpt);
        riadat->abs = (double *) malloc(sizeof(double) * riadat->NBpt);

        for(long i = 0; i < riadat->NBpt; i++)
        {
            int r = fscanf(fp, "%lf %lf %lf\n", &v0, &v1, &v2);
            if(r == 3)
            {
                riadat->lambda[i] = v0;
                riadat->rindex[i] = v1;
                riadat->abs[i] = v2;
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


ATMOSPHERE_MODEL AtmosphereModel_Create_from_CONF(
    const char *restrict CONFFILE,
    float slambda
)
{
    char KEYWORD[200];
    char CONTENT[200];
    
    FILE *fp;

	ATMOSPHERE_MODEL atm;

	DEBUG_TRACEPOINT("Reading atmosphere parameters from file %s", CONFFILE);

    strcpy(KEYWORD, "ZENITH_ANGLE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.ZenithAngle = atof(CONTENT);


    strcpy(KEYWORD, "TIME_DAY_OF_YEAR");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.TimeDayOfYear = atoi(CONTENT);

    strcpy(KEYWORD, "TIME_LOCAL_SOLAR_TIME");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.TimeLocalSolarTime = atoi(CONTENT);


    strcpy(KEYWORD, "SITE_LATITUDE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteLat = atof(CONTENT);

    strcpy(KEYWORD, "SITE_LONGITUDE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteLong = atof(CONTENT);


    strcpy(KEYWORD, "SITE_ALT");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteAlt = atof(CONTENT);



    // temperature at site
    strcpy(KEYWORD, "SITE_TP_AUTO");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteTPauto = atoi(CONTENT);

    // temperature at site
    strcpy(KEYWORD, "SITE_TEMP");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteTemp = atof(CONTENT);

    // pressure at site
    strcpy(KEYWORD, "SITE_PRESS");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SitePress = atof(CONTENT);





    strcpy(KEYWORD, "CO2_PPM");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.CO2_ppm = atof(CONTENT);


    // water
    strcpy(KEYWORD, "SITE_H2O_METHOD");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteH2OMethod = atoi(CONTENT);

    strcpy(KEYWORD, "SITE_TPW");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteTPW = atof(CONTENT);

    strcpy(KEYWORD, "SITE_RH");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SiteRH = atof(CONTENT);

    strcpy(KEYWORD, "SITE_PW_SCALEH");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    atm.SitePWSH = atof(CONTENT);



    DEBUG_TRACEPOINT("Building standard atmoshpere model");
    AtmosphereModel_build_stdAtmModel(atm, "atm.txt");


	
	atm.speciesRIA.NBspecies = atmNBspecies;

	DEBUG_TRACEPOINT("load refractive index data from %d species", atm.speciesRIA.NBspecies);

	atm.speciesRIA.RIA_species = (RIAdata*) malloc(sizeof(RIAdata)*atm.speciesRIA.NBspecies);
	
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
	
	
	for(int sp=0; sp<atm.speciesRIA.NBspecies; sp++)
	{
		char RIAfname[200];
		sprintf(RIAfname, "./RefractiveIndices/RIA_%s.dat", atm.speciesRIA.RIA_species[sp].name);
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

    atm.speciesRIA.RIA_species[speciesH2].Z =  1.0006;
    // 1.013 bar and 15 °C
    // ref: http://encyclopedia.airliquide.com

    atm.speciesRIA.RIA_species[speciesO3].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesN].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesO].Z = 1.0; // TO BE UPDATED

    atm.speciesRIA.RIA_species[speciesH].Z = 1.0; // TO BE UPDATED

	






    // ***************** refractive index as a function of lambda at site ****************************

	DEBUG_TRACEPOINT("Write RindexSite.txt");
    fp  = fopen("RindexSite.txt", "w");
    for(double lambda = 0.2e-6; lambda < 20.0e-6; lambda *= 1.0 + 1e-6)
    {	
		RIAvalue ria = AtmosphereModel_stdAtmModel_ria(atm, atm.SiteAlt, lambda, 0);        
        fprintf(fp, "%.8g %.14f %.14f\n", lambda, ria.rindex, ria.abscoeff);
        //  printf("%g %.10f\n", lambda, n);
    }
    fclose(fp);


	DEBUG_TRACEPOINT("Write Refract.txt");
    fp = fopen("Refract.txt", "w");
    fclose(fp);


    // *************** LIGHT PATH THROUGH ATMOSPHERE AT TWO SEPARATE WAVELENGTHS ***********************






    AtmosphereModel_RefractionPath(atm, 0.55e-6, atm.ZenithAngle, 1);
	EXECUTE_SYSTEM_COMMAND("mv refractpath.txt refractpath_0550.txt");

    AtmosphereModel_RefractionPath(atm, slambda, atm.ZenithAngle, 1);
    EXECUTE_SYSTEM_COMMAND("mv refractpath.txt refractpath_%04ld.txt",
            (long)(1e9 * slambda + 0.5));


    // refractive angle as a function of wavelength
    printf("Computing refraction angle ... \n");
    fflush(stdout);
    fp = fopen("RefractAngle.dat", "w");
    for(double l = 0.5e-6; l < 10.0e-6; l *= 1.0 + 1e-3)
    {
        fprintf(fp, "%20.18f  %.6f\n", l, AtmosphereModel_RefractionPath(atm, l, atm.ZenithAngle,
                0));
    }
    fclose(fp);
    printf("done\n");
    fflush(stdout);

    // **************** REFRACTION AND TRANSMISSION AS A FUNCTION OF WAVELENGTH *********************

    // precompute wavelength array and indices
    /*  llistep = 1;
      llistart = 0;
      while(RIA_N2_lambda[llistart]<4.1e-6)
          llistart++;
      lliend = llistart;
      while(RIA_N2_lambda[lliend]<4.2e-6)
          lliend++;
      NB_comp_array = (lliend-llistart)/llistep;
      comp_array_lambda = (double*) malloc(sizeof(double)*NB_comp_array);
      comp_array_lli = (long*) malloc(sizeof(double)*NB_comp_array);

      for(li=0; li<NB_comp_array; li++)
      {
          lli = llistart + li*llistep;
          comp_array_lambda[li] = RIA_N2_lambda[lli];
          comp_array_lli[li] = lli;
      }

      for(li=0; li<NB_comp_array; li++)
      {
          lambda = comp_array_lambda[li];
          lliprecomp = comp_array_lli[li];
          rangle = AtmosphereModel_RefractionPath(lambda, ZenithAngle, 0);
          fp = fopen("Refract.txt", "a");
          fprintf(fp, "%g %.12f %.12f\n", lambda, rangle, v_TRANSM);
          fclose(fp);
      }
      lliprecomp = -1;

      free(comp_array_lambda);
      free(comp_array_lli);
    */

    return atm;
}
