#include "elp82b.h"

#include <math.h>
#include <stdio.h>

#define epochJ2000(j) (((j) - 2451545.0) / 36525)  /* Julian centuries since J2000 */

#define fixangle(a) ((a) - 360.0 * (floor((a) / 360.0)))
#define fixhour(a) ((a) - 24.0 * (floor((a) / 24.0)))

#define abs(x) (((x) >= 0) ? (x) : (-(x)))
#define abshour(x) ((fixhour(x) <= 12) ? fixhour(x) : 24 - fixhour(x))

#define torad(d) ((d) * (M_PI / 180.0))
#define todeg(d) ((d) * (180.0 / M_PI))

#define dsin(x) (sin(torad((x))))
#define dcos(x) (cos(torad((x))))
#define dtan(x) (tan(torad((x))))

double RMS = 0;
int totaltests = 0;

/***************************************************************/
/*                                                             */
/*  MoonRise --  Calculate time of moonrise/set:               */
/*                                                             */
/*  Algorithm from Explanatory supplement to the astronomical  */
/*  almanac, from The nautical almanac office of the US naval  */
/*  observatory, 1992, pp. 486-489.                            */
/*                                                             */
/***************************************************************/
double MoonRise(int rise, double jd, double lat, double lon, double tz)
{
    double ec, ec0, d0;
    double gst, gha;
    double dec, ra, cosra;
    double mp[3], l, b, e;
    double par, a, ut, ut0, cost;
    int ncosfail, ndayfail;

    ec = epochJ2000(jd);
    e = 23.4392917 - 0.013005 * ec; /* Mean obliquity of the ecliptic. */

    ut = 12.0; /* First iteration */
    ncosfail = 0;
    ndayfail = 0;

    do {
    ut0 = ut;

    ec0 = ec + (ut0 - tz) / 876600.0; /* Julian centuries since J2000 */
    d0 = ec0 * 36525;                 /* Julian days since J2000 */

    /* Calculate the moon's position */
    GetElp82bSphericalCoor(ec0, mp);
    l = mp[0] / 3600.0;
    b = mp[1] / 3600.0;

    /* Convert from ecliptic to equatorial co-ordinates */
    dec = todeg(asin( dsin(b) * dcos(e) + dcos(b) * dsin(e) * dsin(l) ));
    cosra = dcos(l) * dcos(b) / dcos(dec);

    /* Check the sign of sin(ra), to determine which branch of acos to use */
    if (dcos(b) * dsin(l) * dcos(e) >= dsin(b) * dsin(e))
        ra = todeg(acos(cosra));
    else
        ra = -todeg(acos(cosra));

    /* Greenwich mean sidereal time in hours */
    gst = 18.697374558 + 24.06570982441908 * d0;
    /* Greenwich hour angle of the moon in hours */
    gha = gst - ra / 15.0;

    /* Horizontal parallax of the moon.*/
    par = todeg(atan( 6378.137 / mp[2] ));

    /* Apparent altitude of upper limb of moon at moonrise/set */
    a = -34.0 / 60.0 + 0.7275 * par;

    /* Assuming the moon remains stationary, get the local hour angle of rise/set */
    cost = (dsin(a) - dsin(lat) * dsin(dec)) / (dcos(lat) * dcos(dec));

    /* Check if the moon never rises/sets */
    if (cost > 1) {
        cost = 1;  /* Never sets */
        if (ncosfail < 0) ncosfail = 0;
        ncosfail++;
    } else if (cost < -1) {
        cost = -1;  /* Never rises */
        if (ncosfail > 0) ncosfail = 0;
        ncosfail--;
    } else
        ncosfail = 0;

    /* Calculate UTC time of moonrise/set. */
    if (rise == 1)
        ut = ut0 - gha - ( lon + todeg(acos(cost)) ) / 15;
    else
        ut = ut0 - gha - ( lon - todeg(acos(cost)) ) / 15;
    ut = fixhour(ut);

    /* Check if we've rolled over to the next/previous day */
    if (abs(ut - ut0) > 12)
        ndayfail++;

    /* Debugging
    printf("Lat %.2f, Lon %.2f, Days %.1f, ut0 %09f, GHA %09f, t %f, ut %09f\n",
            lat, lon, d0, ut0, fixhour(gha), todeg(acos(cost)), ut);
    */

    /* Fail with no rise/set */
    if (abs(ncosfail) > 3 || ndayfail > 1)
        return -1;

    } while (abs(ut - ut0) >= 0.008);

    return ut;
}

void hourToClock(double t, char* result)
{
    int hour, min;

    hour = floor(t);
    min = floor((t - hour) * 60 + 0.5);

    sprintf(result, "%02dh%02d", hour, min);
}

double clockToHour(double t)
{
    double hour, min;

    hour = floor(t / 100);
    min = floor(t - hour * 100);

    return hour + min / 60.0;
}

int test(int rise, double jd, char* date, double lat, double lon,
         double tz, double target, int show)
{
    double result, t, err;
    char sresult[5], starget[5];

    double maxerr = 0.08;

    result = MoonRise(rise, jd, lat, lon, tz);

    if (target >= 0) {
        t = clockToHour(target);
        hourToClock(t, starget);
    } else {
        t = -1;
        sprintf(starget, "None");
    }

    if (result >= 0)
        hourToClock(result, sresult);
    else
        sprintf(sresult, "None");

    err = abshour(result - t);
    RMS += err * err;
    totaltests++;

    if (err > maxerr || show)
    {
        printf("%s at %s (moon%s), Lat %02.2f, Lon %03.2f, TZ %+.1f\n",
               err > maxerr ? "Fail" : "Success", date, rise ? "rise": "set",
               lat, lon, tz);
        printf("  Result %s, Target %s\n", sresult, starget);
        /*printf("  Result %f, Target %f\n", result, t);*/
    }

    return err <= maxerr;
}

int testsuite(int show)
{
    double lat, lon, tz;

    RMS = 0;
    totaltests = 0;

    lat = 0; lon = 0; tz = 0;

    test(1, 2451547.5, "2000 Jan 04", lat, lon, tz,  405, show);
    test(0, 2451547.5, "2000 Jan 04", lat, lon, tz, 1628, show);
    test(1, 2451548.5, "2000 Jan 05", lat, lon, tz,  452, show);
    test(0, 2451548.5, "2000 Jan 05", lat, lon, tz, 1716, show);
    test(1, 2451549.5, "2000 Jan 06", lat, lon, tz,  541, show);
    test(0, 2451549.5, "2000 Jan 06", lat, lon, tz, 1804, show);

    lat =  -35.0 - 18.0 / 60.0; /* - 27.0 / 3600.0; Latitude of observer */
    lon =  149.0 + 07.0 / 60.0; /* + 27.9 / 3600.0; Longitude of observer */
    tz = 11;

    test(1, 2457771.5, "2017 Jan 18", lat, lon, tz, 2359, show);
    test(0, 2457771.5, "2017 Jan 18", lat, lon, tz, 1143, show);
    test(1, 2457772.5, "2017 Jan 19", lat, lon, tz,   -1, show);
    test(0, 2457772.5, "2017 Jan 19", lat, lon, tz, 1239, show);
    test(1, 2457773.5, "2017 Jan 20", lat, lon, tz,   31, show);
    test(0, 2457773.5, "2017 Jan 20", lat, lon, tz, 1334, show);

    printf("RMS: %.2f minutes\n", sqrt(RMS / totaltests) * 60.0);

    return 0;
}

int main()
{
    /*
    double mp[3];
    GetElp82bSphericalCoor(0, mp);
    printf("0: %f, 1: %f, 2: %f\n", mp[0] / 3600.0, mp[1] / 3600.0, mp[2]);
    */
    return testsuite(1);
}
