/**
 * @file    AtmosphereModel_stdAtmModel_ria.h
 *
 */

#ifndef _ATMOSPHEREMODEL_STDATMMODEL_RIA_H
#define _ATMOSPHEREMODEL_STDATMMODEL_RIA_H

RIAvalue AtmosphereModel_stdAtmModel_ria(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
);

double AtmosphereModel_stdAtmModel_N(
    ATMOSPHERE_MODEL atm,
    float alt,
    float lambda,
    int mode
);

#endif
