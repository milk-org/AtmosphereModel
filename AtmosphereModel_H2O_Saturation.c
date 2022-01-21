/**
 * @file    AtmosphereModel_H2O_Saturation.c
 *
 */

#include <math.h>

/// using " The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use ", Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31, Issue 2, pp. 387535
/// output : saturation H2O pressure [Pa]

double AtmosphereModel_H2O_Saturation(double T)
{
    double val;
    double ups;

    double Tc = 647.096;    // critical temperature [K]
    double Pc = 22064000.0; // critical pressure [Pa]
    //   double P0 = 101325.0; // standard pressure [Pa]

    double C1 = -7.85951783;
    double C2 = 1.84408259;
    double C3 = -11.7866497;
    double C4 = 22.6807411;
    double C5 = -15.9618719;
    double C6 = 1.80122502;

    double Tn = 273.16; // tripple point
    double a0 = -13.928169;
    double a1 = 34.707823;
    double Pn = 611.657;

    if (T > 273.0)
    {
        ups = 1.0 - T / Tc;
        val = Pc * exp(Tc / T *
                       (C1 * ups + C2 * pow(ups, 1.5) + C3 * pow(ups, 3.0) +
                        C4 * pow(ups, 3.5) + C5 * pow(ups, 4.0) +
                        C6 * pow(ups, 7.5)));
    }
    else
    {
        ups = T / Tn;
        val = Pn *
              exp(a0 * (1.0 - pow(ups, -1.5)) + a1 * (1.0 - pow(ups, -1.25)));
    }

    return val;
}
