/**
 * @file    AtmosphereModel_build_stdAtmModel.c
 *
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "AtmosphereModel.h"
#include "AtmosphereModel_H2O_Saturation.h"
#include "AtmosphereModel_save_stdAtmModel.h"
#include "nrlmsise-00.h"


int AtmosphereModel_build_stdAtmModel(
    ATMOSPHERE_MODEL atm,
    const char *restrict fname
)
{
    FILE *fp;
    struct nrlmsise_output output[10000]; // 100m steps from h=0 to h=100km, 10m steps
    struct nrlmsise_input input[10000];
    struct nrlmsise_flags flags;
    //	struct ap_array aph; // magnetic values
    int i;
    //int j;
    double TotPart; // total number of particules per cm2
    double TotPart1, TotPart2; // for H2O
    double h; // elevation
    double X;
    double coeff;

    //double LLN2, LLO2, LLCO2, LLHe, LLAr, LLH;

    double LoschmidtConstant =  2.6867805e25;
    double *TotPart0;
    double densH20_site;
    double H2OSatTempSite;
    long i0;
    float ifrac;

    double tpw, tpw0, pwsh, deltapwsh;
    long iter;
    int dir, odir;


    int k;
    double deltah = 0.0;
    double PressCoeff = 1.0;


    /* input values */
    //  	for (i=0;i<7;i++)
    //		aph.a[i]=100;
    flags.switches[0] = 0;
    for(i = 1; i < 24; i++)
    {
        flags.switches[i] = 1;
    }

    for(i = 0; i < 10000; i++)
    {
        input[i].doy = atm.TimeDayOfYear;
        input[i].year = 0; /* without effect */
        input[i].alt = 0.01 * i; // [km]
        input[i].g_lat = atm.SiteLat;
        input[i].g_long = atm.SiteLong;
        input[i].lst = atm.TimeLocalSolarTime;
        input[i].sec = (int)(3600.0 * (atm.TimeLocalSolarTime - atm.SiteLong / 15.0));
        if(input[i].sec < 0)
        {
            input[i].sec += 3600 * 24;
        }
        input[i].f107A = 150;
        input[i].f107 = 150;
        input[i].ap = 4;
    }
    fp = fopen(fname, "w");
    /* evaluate 0 to 1000 */
    fprintf(fp, "#  alt[m]  N2  O2  Ar  H  He  O   N  density   Temperature\n");
    TotPart = 0.0;
    TotPart1 = 0.0;
    TotPart2 = 0.0;
    //dens0 = 1;

    i = 0;
    printf("Site altitude: %f m\n", atm.SiteAlt);



    h = atm.SiteAlt;
    input[i].alt = 0.001 * h;
    gtd7(&input[i], &flags, &output[i]);
    atm.vprof_dens_species[speciesNe].val[i] = output[i].d[2] * 2.328e-5;
    
    atm.vprof_denstot.val[i] = output[i].d[0] + output[i].d[1] + output[i].d[2] + output[i].d[3] +
                 output[i].d[4] + output[i].d[6] + output[i].d[7] + output[i].d[8] + atm.vprof_dens_species[speciesNe].val[i];
    
    atm.vprof_pressure.val[i] = atm.vprof_denstot.val[i] * 1.0e6 / LoschmidtConstant * (output[i].t[1] /
                  273.15);
    printf("Temperature at site altitude (%12f m) = %8f K\n", atm.SiteAlt,
           output[i].t[1]);
    printf("Pressure at site altitude    (%12f m) = %8f atm\n", atm.SiteAlt,
           atm.vprof_pressure.val[i]);



    if(atm.SiteTPauto == 0) // adjust temperature and pressure
    {
        // search for altitude to match site temperature
        h = atm.SiteAlt;
        for(k = 0; k < 10; k++)
        {
            input[i].alt = 0.001 * h;
            gtd7(&input[i], &flags, &output[i]);
            atm.vprof_dens_species[speciesNe].val[i] = output[i].d[2] * 2.328e-5;
            atm.vprof_denstot.val[i] = output[i].d[0] + output[i].d[1] + output[i].d[2] + output[i].d[3] +
                         output[i].d[4] + output[i].d[6] + output[i].d[7] + output[i].d[8] + atm.vprof_dens_species[speciesNe].val[i];
            atm.vprof_pressure.val[i] = atm.vprof_denstot.val[i] * 1.0e6 / LoschmidtConstant * (output[i].t[1] /
                          273.15);

            printf("h = %.0f m   Temp = %f  (%f)  Press = %f  (%f)\n", h, output[i].t[1],
                   atm.SiteTemp, atm.vprof_pressure.val[i], atm.SitePress);
            h += 1000.0 * (output[i].t[1] - atm.SiteTemp) / 6.5;
        }
        deltah = h - atm.SiteAlt;
        PressCoeff = atm.SitePress / atm.vprof_pressure.val[i];
        printf("deltah = %f m   PressCoeff = %f\n", deltah, PressCoeff);
    }
    else
    {
        deltah = 0.0;
        PressCoeff = 1.0;
    }




    TotPart0 = (double *) malloc(sizeof(double) *
                                 10000); // total number of particles per cm3

    for(i = 0; i < 10000; i++)
    {
        h = 10.0 * i + deltah; // 10 m step
        input[i].alt = 0.001 * h;
        gtd7(&input[i], &flags, &output[i]);
        fprintf(fp, "%6.0f %12f %12f %12f %12f %12f %12f %12f %.8g %5f\n",
                input[i].alt * 1000.0, output[i].d[2], output[i].d[3], output[i].d[4],
                output[i].d[6], output[i].d[0], output[i].d[1] + output[i].d[8], output[i].d[7],
                output[i].d[5], output[i].t[1]);

        output[i].d[0] *= PressCoeff;
        output[i].d[1] *= PressCoeff;
        output[i].d[2] *= PressCoeff;
        output[i].d[3] *= PressCoeff;
        output[i].d[4] *= PressCoeff;
        output[i].d[5] *= PressCoeff;
        output[i].d[6] *= PressCoeff;
        output[i].d[7] *= PressCoeff;
        output[i].d[8] *= PressCoeff;

        TotPart0[i] = 0.0;
        TotPart0[i] += output[i].d[0];
        TotPart0[i] += output[i].d[1];
        TotPart0[i] += output[i].d[2];
        TotPart0[i] += output[i].d[3];
        TotPart0[i] += output[i].d[4];
        TotPart0[i] += output[i].d[6];
        TotPart0[i] += output[i].d[7];
        TotPart0[i] += output[i].d[8];


        atm.vprof_dens_species[speciesNe].val[i] = output[i].d[2] * 2.33e-5;
        TotPart0[i] += atm.vprof_dens_species[speciesNe].val[i]; // total density

        atm.vprof_dens_species[speciesCH4].val[i] = TotPart0[i] * 2e-6;
        if(h < 45000.0)
        {
            atm.vprof_dens_species[speciesCH4].val[i] *= 1.0 - h / 45000.0;
        }
        else
        {
            atm.vprof_dens_species[speciesCH4].val[i] = 0.0;
        }
        TotPart0[i] += atm.vprof_dens_species[speciesCH4].val[i]; // total density

        atm.vprof_dens_species[speciesKr].val[i] = output[i].d[2] * 1.46e-6;
        TotPart0[i] += atm.vprof_dens_species[speciesKr].val[i]; // total density

        atm.vprof_dens_species[speciesH2].val[i] = output[i].d[2] * 7.04e-7;
        TotPart0[i] += atm.vprof_dens_species[speciesH2].val[i]; // total density

        atm.vprof_dens_species[speciesO3].val[i] = 300.0 * 2.69e20 / (4250.0 * sqrt(2.0 * M_PI)) * exp(-0.5 * pow(((
                        h - 25000.0) / 4250.0), 2.0)); // [m-3]
        atm.vprof_dens_species[speciesO3].val[i] *= 1.0e-6;



        if(h < 70000)
        {
            atm.vprof_dens_species[speciesCO2].val[i] = atm.CO2_ppm * 1e-6 * TotPart0[i];
        }
        else
        {
            atm.vprof_dens_species[speciesCO2].val[i] = (atm.CO2_ppm - 0.007 * (h - 70000)) * 1e-6 * TotPart0[i];
        }

        output[i].d[3] -= atm.vprof_dens_species[speciesCO2].val[i];


        //if(i == 0)
        //{
        //    dens0 = TotPart0[i];
        //}
        TotPart0[i] *= 1000.0; // convert from cm-3 to cm-2   each bin is 10m = 1000cm
        // TotPart0 is the number of particles per cm2 for the altitude bin

        atm.vprof_temperature.val[i] = output[i].t[1];

        if(h > atm.SiteAlt)
        {
            TotPart += TotPart0[i]; // cumulative number of particles per cm2
            TotPart1 += 2.5e-6 * TotPart0[i] * exp(-3.0 * pow(h / 100000.0,
                                                   4.0)); // fixed 2.5e-6 mixing ratio term
            TotPart2 += (exp(-h / atm.SitePWSH)) * TotPart0[i] * exp(-3.0 * pow(h / 100000.0,
                        4.0)); // exponential mixing term (not yet scaled)
            //			printf("===== %f %f %f\n", TotPart, TotPart0, TotPart1);
        }

        // for now, assume fixed background mixing ratio
        atm.vprof_dens_species[speciesH2O].val[i] = 2.5e-6 * TotPart0[i] / 1000.0; // cm-3
    }
    fclose(fp);

    printf("TOTAL particles : %g part.cm^-2\n", TotPart);
    printf("TOTAL particles H2O fixed mixing ratio : %g part.cm^-2    ( SitePWSH = %f, SiteAlt = %f m)\n",
           TotPart1, atm.SitePWSH, atm.SiteAlt);


    if(atm.SiteH2OMethod ==
            1) // TPW, SH -> adjust exponential coeff to meet total precipitable water
    {
        X = (atm.SiteTPW * 0.1) * 6.0221413e23 /
            18.0; // total particules of H2O per cm2 = water mass [g] x avogadro's nuber / 18.0
        printf("TOTAL H2O particles: %g part.cm^-2\n", X);

        atm.alpha1H2O = (X - TotPart1) / TotPart2;
        if(atm.alpha1H2O < 0.0)
        {
            printf("ERROR: total precititable water value is too low, and not consistent with atmospheric model\n");
            exit(0);
        }
        printf("alpha1H2O = %g\n", atm.alpha1H2O);

        for(i = 0; i < 10000; i++)
        {
            h = 10.0 * i;
            atm.vprof_dens_species[speciesH2O].val[i] += atm.alpha1H2O * (exp(-h / atm.SitePWSH)) * (TotPart0[i] / 1000.0) * exp(
                              -3.0 * pow(h / 100000.0, 4.0));
        }
    }

    if(atm.SiteH2OMethod == 2) // RH, SH -> compute total precipitable water
    {
        i0 = (long)(atm.SiteAlt / 10.0);
        ifrac = atm.SiteAlt / 10.0 - i0;

        H2OSatTempSite = AtmosphereModel_H2O_Saturation((1.0 - ifrac) * atm.vprof_temperature.val[i0]
                         + ifrac * atm.vprof_temperature.val[i0 + 1]); // [Pa]

        densH20_site = (atm.SiteRH / 100.0) * LoschmidtConstant * 1e-6 * H2OSatTempSite /
                       101325.0; // cm3
        //	RH[i] = (densH2O[i]*1e6/LoschmidtConstant)*101325.0/AtmosphericTurbulence_H2O_Saturation(temperature[i]);
        TotPart2 = 0.0;
        for(i = 0; i < 10000; i++)
        {
            h = 10.0 * i;
            //            densH2O[i] /= 1000.0;
            atm.vprof_dens_species[speciesH2O].val[i] += densH20_site * exp(-(h - atm.SiteAlt) / atm.SitePWSH) *
                          (TotPart0[i] / TotPart0[i0]) * exp(-3.0 * pow(h / 100000.0, 4.0));
            if(h > atm.SiteAlt)
            {
                TotPart2 += 100.0 * densH20_site * exp(-(h - atm.SiteAlt) / atm.SitePWSH) * exp(
                                -3.0 * pow(h / 100000.0, 4.0));
            }
        }

        printf("Total precipitable water : %g mm\n",
               (TotPart1 + TotPart2) * 10.0 * 18.0 / 6.022e23);
    }


    if(atm.SiteH2OMethod == 3) // RH, TPW -> compute SH
    {

        X = (atm.SiteTPW * 0.1) * 6.0221413e23 /
            18.0; // total particules of H2O per cm2 = water mass [g] x avogadro's nuber / 18.0
        printf("TOTAL H2O particles: %g part.cm^-2\n", X);

        i0 = (long)(atm.SiteAlt / 10.0);
        ifrac = atm.SiteAlt / 10.0 - i0;

        H2OSatTempSite = AtmosphereModel_H2O_Saturation((1.0 - ifrac) * atm.vprof_temperature.val[i0]
                         + ifrac * atm.vprof_temperature.val[i0 + 1]); // [Pa]
        densH20_site = (atm.SiteRH / 100.0) * LoschmidtConstant * 1e-6 * H2OSatTempSite /
                       101325.0; // cm3

        tpw = 0.0;
        for(i = 0; i < 10000; i++)
        {
            h = 10.0 * i;
            if(h > atm.SiteAlt)
            {
                tpw += atm.vprof_dens_species[speciesH2O].val[i] * 1000.0; // particules per cm2
            }
        }
        printf("    H2O :   %g particle/cm2\n", tpw);
        tpw /= 6.0221413e23;
        tpw *= 18.0;
        tpw *= 10.0;
        printf("    H2O :   %f mm TPW\n", tpw);


        printf("missing : %f mm\n", atm.SiteTPW - tpw);
        tpw0 = atm.SiteTPW - tpw;

        pwsh = 1000.0;
        deltapwsh = 100.0;
        dir = 1;
        for(iter = 0; iter < 15; iter++)
        {
            TotPart2 = 0.0;
            for(i = 0; i < 10000; i++)
            {
                h = 10.0 * i;
                if(h > atm.SiteAlt)
                {
                    TotPart2 += 1000.0 * densH20_site * exp(-(h - atm.SiteAlt) / pwsh) * exp(-3.0 * pow(
                                    h / 100000.0, 4.0));
                }
            }
            tpw = TotPart2 * 10.0 * 18.0 / 6.022e23;
            printf("SH = %8f m   -> Total precipitable water term 2: %g mm\n", pwsh, tpw);

            printf("  %f  %f %f  %f\n", pwsh, tpw, tpw0, deltapwsh);

            odir = dir;
            if(tpw > tpw0) // too much water
            {
                pwsh -= deltapwsh;
                dir = -1;
            }
            else // not enough
            {
                pwsh += deltapwsh;
                dir = 1;
            }

            if(odir != dir)
            {
                deltapwsh *= 0.5;
            }
        }

        for(i = 0; i < 10000; i++)
        {
            h = 10.0 * i;
            //          densH2O[i] /= 1000.0;
            atm.vprof_dens_species[speciesH2O].val[i] += densH20_site * exp(-(h - atm.SiteAlt) / pwsh) *
                          (TotPart0[i] / TotPart0[i0]) * exp(-3.0 * pow(h / 100000.0, 4.0));
        }

    }


    // Checking H2O
    // measuring TPW
    tpw = 0.0;
    for(i = 0; i < 10000; i++)
    {
        h = 10.0 * i;
        if(h > atm.SiteAlt)
        {
            tpw += atm.vprof_dens_species[speciesH2O].val[i] * 1000.0; // particules per cm2
        }
    }
    printf("    H2O :   %g particle/cm2\n", tpw);
    tpw /= 6.0221413e23;
    tpw *= 18.0;
    tpw *= 10.0;
    printf("    H2O :   %f mm TPW\n", tpw);

    i0 = (long)(atm.SiteAlt / 10.0);
    ifrac = atm.SiteAlt / 10.0 - i0;

    H2OSatTempSite = AtmosphereModel_H2O_Saturation((1.0 - ifrac) * atm.vprof_temperature.val[i0]
                     + ifrac * atm.vprof_temperature.val[i0 + 1]);
    printf("    RH  :   %f %%\n",
           atm.vprof_dens_species[speciesH2O].val[i0] / LoschmidtConstant * 1e6 / H2OSatTempSite * 101325.0 * 100.0);




    //fp = fopen(fname, "w");
    /* evaluate 0 to 1000 */
    //  fprintf(fp, "#  1:alt[m]  2:denstot[part/cm3] 3:N2  4:O2  5:Ar  6:H2O  7:CO2  8:Ne  9:He  10:CH4  11:Kr  12:H2  13:N  14:O  15:H    16:density[g/cm3]  17:temperature[K] 18:pressure[stdatm]  19:RH\n");
    TotPart = 0.0;
    TotPart1 = 0.0;
    TotPart2 = 0.0;
    for(i = 0; i < 10000; i++)
    {
        h = 10.0 * i;
        //	gtd7(&input[i], &flags, &output[i]);

        if(h > atm.SiteAlt)
        {
            TotPart1 += atm.vprof_dens_species[speciesH2O].val[i] * 1000.0;
        }

        // total particule densisty without H2O and O3
        atm.vprof_denstot.val[i] = output[i].d[0] + output[i].d[1] + output[i].d[2] + output[i].d[3] +
                     output[i].d[4] + output[i].d[6] + output[i].d[7] + output[i].d[8] + atm.vprof_dens_species[speciesCO2].val[i] +
                     atm.vprof_dens_species[speciesNe].val[i];

        coeff = (atm.vprof_denstot.val[i] + atm.vprof_dens_species[speciesH2O].val[i] + atm.vprof_dens_species[speciesO3].val[i]) / atm.vprof_denstot.val[i];
        //			(output[i].d[0]+output[i].d[1]+output[i].d[2]+output[i].d[3]+output[i].d[4]+output[i].d[6]+output[i].d[7]+output[i].d[8]+densCO2[i]+densNe[i]+densH2O[i])/(output[i].d[0]+output[i].d[1]+output[i].d[2]+output[i].d[3]+output[i].d[4]+output[i].d[6]+output[i].d[7]+output[i].d[8]+densCO2[i]+densNe[i]);

        //printf("%10f  H2O coeff = %f\n", h, coeff);

        output[i].d[0] /= coeff;  	// He
        atm.vprof_dens_species[speciesHe].val[i] = output[i].d[0];

        output[i].d[1] /= coeff;	// O
        atm.vprof_dens_species[speciesO].val[i] = output[i].d[1];

        output[i].d[2] /= coeff;	// N2
        atm.vprof_dens_species[speciesN2].val[i] = output[i].d[2];

        output[i].d[3] /= coeff;	// O2
        atm.vprof_dens_species[speciesO2].val[i] = output[i].d[3];

        output[i].d[4] /= coeff;	// Ar
        atm.vprof_dens_species[speciesAr].val[i] = output[i].d[4];

        output[i].d[6] /= coeff;	// H
        atm.vprof_dens_species[speciesH].val[i] = output[i].d[6];

        output[i].d[7] /= coeff;	// N
        atm.vprof_dens_species[speciesN].val[i] = output[i].d[7];

        output[i].d[8] /= coeff;	// anomalous O
        atm.vprof_dens_species[speciesO].val[i] += output[i].d[8];

        atm.vprof_dens_species[speciesH2O].val[i] /= coeff;

        atm.vprof_dens_species[speciesCO2].val[i] /= coeff;

        atm.vprof_dens_species[speciesNe].val[i] /= coeff;

        atm.vprof_dens_species[speciesKr].val[i] = atm.vprof_dens_species[speciesN2].val[i] / 78.084 * 1.14e-6;
        atm.vprof_dens_species[speciesH2].val[i] = atm.vprof_dens_species[speciesN2].val[i] / 78.084 * 5.5e-7;

        atm.vprof_density.val[i] = output[i].d[5];
        atm.vprof_temperature.val[i] = output[i].t[1];
        atm.vprof_pressure.val[i] = atm.vprof_denstot.val[i] * 1.0e6 / LoschmidtConstant * (output[i].t[1] /
                      273.15);

        atm.vprof_RH.val[i] = (atm.vprof_dens_species[speciesH2O].val[i] * 1e6 / LoschmidtConstant) * 101325.0 /
                AtmosphereModel_H2O_Saturation(atm.vprof_temperature.val[i]);


        //        fprintf(fp, "%6.0f   %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g     %.8g    %8.3lf   %12.10lf  %.5f\n", input[i].alt*1000.0, denstot[i], densN2[i], densO2[i], densAr[i], densH2O[i], densCO2[i], densNe[i], densHe[i], densCH4[i], densKr[i], densH2[i], densN[i], densO[i], densH[i],  output[i].d[5], temperature[i], pressure[i], RH[i]);
    }
    //  fclose(fp);
    printf("TOTAL = %g\n", TotPart1);

    free(TotPart0);


    AtmosphereModel_save_stdAtmModel(atm, fname);

    return(0);
}
