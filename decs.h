#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define FOV 128.
#define NX  256
#define NY  256

#define RAINBOW 0
#define AFMHOT 1
#define BW 0

#define NDIM  4

/* imaging */
void make_ppm(double p[NX][NY], double freq, char filename[]) ;
void rainbow_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
void monika_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
void afmhot_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;

/*find index on a map*/
void xtoij(double x, double y, int *i, int *j, double *delx, double *dely);
/*interpolate vector component or scalar to x,y on the map*/
double interp_scalar2d(double x, double y, double var[NX][NY]);
double hanning_windowed(double a,double b);



