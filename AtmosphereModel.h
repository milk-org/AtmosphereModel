/**
 * @file    AtmosphereModel.h
 * @brief   Function prototypes for AtmosphereModel functions
 *
 * Function prototypes are added here for use by other modules.\n
 * Do not use this .h file for forward declaration internal to the module.\n
 */



#ifndef _ATMOSPHEREMODEL_H
#define _ATMOSPHEREMODEL_H


typdedef struct
{
    float ZenithAngle;
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
} ATMOSPHERE_MODEL;



//#include "create_example_image.h"
//#include "stream_process_loop_simple.h"

#endif
