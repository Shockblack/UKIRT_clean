#include <math.h>

// Function to integrate in scipy's LowLevelCallable format

// c[0] = EWRC; c[1] = B; c[2] = M_RC; c[3] = sigma_RC

double f(int n, double *x, void *user_data) {

    double EWRC = ((double *)user_data)[0];
    double B = ((double *)user_data)[1];
    double M_RC = ((double *)user_data)[2];
    double sigma_RC = ((double *)user_data)[3];

    double EWRGBB = 0.201*EWRC;
    double M_RGBB = M_RC + 0.737;
    double sigma_RGBB = sigma_RC;

    double x1 = exp(B*(x[0]-M_RC));
    double x2 = EWRC * exp(-pow((x[0]-M_RC),2)/(2*pow(sigma_RC,2))) / (sigma_RC*sqrt(2*M_PI));
    double x3 = EWRGBB * exp(-pow((x[0]-M_RGBB),2)/(2*pow(sigma_RGBB,2))) / (sigma_RGBB*sqrt(2*M_PI));

    double lum_func = x1 + x2 + x3;

    return lum_func;
}