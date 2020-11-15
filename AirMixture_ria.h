/**
 * @file    AirMixture_ria.h
 *
 */

#ifndef _AIRMIXTURE_RIA_H
#define _AIRMIXTURE_RIA_H

RIAvalue AirMixture_ria(
	ATM_SPECIES_RIADATA speciesRIA,
    double lambda,
    double *densarray
);

#endif
