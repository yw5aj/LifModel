#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

double gaussian_noise(const double mean, const double std)
{
    static int have_spare = 0;
    static double u1, u2, z1, z2;
    if(have_spare)
    {
        have_spare = 0;
        z2 = sqrt(-2. * log(u1)) * sin(2. * M_PI * u2);
        return mean + std * z2;
    }
    have_spare = 1;
    u1 = ((double) (rand() + 1) / (RAND_MAX + 1));
    u2 = ((double) (rand() + 1) / (RAND_MAX + 1));
    z1 = sqrt(-2. * log(u1)) * cos(2. * M_PI * u2);
    return mean + std * z1;
}