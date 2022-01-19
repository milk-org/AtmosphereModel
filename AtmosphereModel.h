/**
 * @file    AtmosphereModel.h
 * @brief   Function prototypes for AtmosphereModel functions
 *
 * Function prototypes are added here for use by other modules.\n
 * Do not use this .h file for forward declaration internal to the module.\n
 */

#ifndef _ATMOSPHEREMODEL_H
#define _ATMOSPHEREMODEL_H

// number of elevation steps in vertical profiles
#define ATM_VPROF_NBSTEP 10000

// vertical step size [m]
#define ATM_VPROF_STEPSIZE 10.0

typedef struct
{
    double rindex;   // refractive index
    double abscoeff; // absorption coefficient [m-1]
} RIAvalue;

// refractive index and absorption data for a signle molecule
typedef struct
{
    char name[8];
    int init;
    double lambdamin; // [m]
    double lambdamax; // [m]
    long NBpt;
    double *lambda;
    double *rindex;
    double *abs;

    double Z; // compressibility factor

} RIAdata;

// vertical profile
typedef struct
{
    //	long  NBpt; // number of smamples
    //	float dz;   // elevation step
    //	int   init;
    double val[ATM_VPROF_NBSTEP];
} ATMvPROF;

// assign index to each species
#define atmNBspecies 14
#define speciesN2 0
#define speciesO2 1
#define speciesAr 2
#define speciesH2O 3
#define speciesCO2 4
#define speciesNe 5
#define speciesHe 6
#define speciesCH4 7
#define speciesKr 8
#define speciesH2 9
#define speciesO3 10
#define speciesN 11
#define speciesO 12
#define speciesH 13

typedef struct
{
    int NBspecies;
    RIAdata *RIA_species; // refractive index and abs coeff
} ATM_SPECIES_RIADATA;

// Atmosphere model at site
//
typedef struct
{
    double ZenithAngle;

    // refractive indices and absorption coefficients
    ATM_SPECIES_RIADATA speciesRIA;

    int TimeDayOfYear;
    float TimeLocalSolarTime;

    float SiteLat;
    float SiteLong;
    float SiteAlt;
    float CO2_ppm;

    int SiteTPauto;
    float SiteTemp;
    float SitePress;

    float SiteH2OMethod;
    float SiteTPW;
    float SiteRH;
    float SitePWSH;
    float alpha1H2O;

    // vertical profiles
    //long  vprofNBpt; // number of smamples
    //float vprofdz;   // elevation step

    // composition (densities)
    ATMvPROF vprof_denstot; // total density
    ATMvPROF vprof_dens_species[atmNBspecies];

    // physical parameters
    ATMvPROF vprof_density;
    ATMvPROF vprof_temperature;
    ATMvPROF vprof_pressure;
    ATMvPROF vprof_RH;

    long NB_comp_array;
    double *comp_array_lambda;
    long *comp_array_lli;
    int lliprecomp; // if >=0, use this index in array

    int lliprecompN2;  // if >=0, use this index in array
    int lliprecompO2;  // if >=0, use this index in array
    int lliprecompAr;  // if >=0, use this index in array
    int lliprecompH2O; // if >=0, use this index in array
    int lliprecompCO2; // if >=0, use this index in array
    int lliprecompNe;  // if >=0, use this index in array
    int lliprecompHe;  // if >=0, use this index in array
    int lliprecompCH4; // if >=0, use this index in array
    int lliprecompKr;  // if >=0, use this index in array
    int lliprecompH2;  // if >=0, use this index in array
    int lliprecompO3;  // if >=0, use this index in array
    int lliprecompO;   // if >=0, use this index in array
    int lliprecompN;   // if >=0, use this index in array
    int lliprecompH;   // if >=0, use this index in array

} ATMOSPHERE_MODEL;

//#include "create_example_image.h"
//#include "stream_process_loop_simple.h"

#endif
