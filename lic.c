#include "decs.h"
#include <omp.h>
#ifdef _OPENMP
#include <omp.h>
#endif


double imageS[NX][NY][NDIM];
double image[NX][NY];
double imageb[NX][NY];
double imageI[NX][NY];
double vx[NX][NY];
double vy[NX][NY];

//xy in uas, pixel center
double x[NX];
double y[NY];

double dx=FOV/NX;
double dy=FOV/NY;
double startx=-FOV/2.;
double stopx=FOV/2.;
double starty=-FOV/2.;
double stopy=FOV/2.;

void main(){

    FILE *fp;
    int i, j;    
    char* fname = "ipole.dat";
    double dum_dbl;
    int dum_int;
    double freq=230.e9;
    
    fprintf(stdout,"%s\n", fname);
    fp = fopen(fname, "r");
    if (fp == NULL) {
    	fprintf(stderr, "unable to open %s\n", fname);
    	exit(1);
    }
    
    for (i = 0; i < NX; i++) {
    	for (j = 0; j < NY; j++) {
    	    fscanf(fp, "%d %d %lf %lf %lf %lf %lf %lf\n",
	    	    &dum_int, &dum_int, &dum_dbl,&imageS[i][j][0],&imageS[i][j][1],&imageS[i][j][2],&imageS[i][j][3],&dum_dbl);
//	    fscanf(fp, "%d %d %lf %lf %lf %lf %lf \n",
//		   &dum_int, &dum_int, &dum_dbl,&imageS[i][j][0],&imageS[i][j][1],&imageS[i][j][2],&imageS[i][j][3]);
    	}
    	fscanf(fp,"\n");
    }
    fclose(fp);
    
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
	    image[i][j] = pow(imageS[i][j][0],0.2);
    make_ppm(image, freq, "ipole_fnu.ppm");

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            image[i][j] = log(imageS[i][j][0] + 1.e-50);
    make_ppm(image, freq, "ipole_lfnu.ppm");

    //in [uas]
    for (int i = 0; i < NX; i++){
	x[i]=-FOV/2. + (i+0.5)*dx;
	y[i]=-FOV/2. + (i+0.5)*dy;	
    }
    
    double amp,ang;
    //construct a vector field for each pixel
    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){
	    amp=1.;//FOV/NX;
	    ang=0.5*atan2(imageS[i][j][2],imageS[i][j][1])+M_PI/2.;
//	    ang=0.5*atan2(0.0,-fabs(imageS[i][j][1]));
	    vx[i][j] = (-sin(ang)*amp);
	    vy[i][j] = ( cos(ang)*amp);
//	    imageI[i][j]=rand();
//	    imageI[i][j]=rand()*imageS[i][j][0];
//	    imageI[i][j]=rand()*log(imageS[i][j][0]+1e-50);
	    imageI[i][j]=rand()*pow(imageS[i][j][0],0.1);
	    image[i][j]=0.0;
	}
    }
    

    // define length scale on the files L in uas
    // we do not integrate further
    double L=3;//10.;//80.; //uas//we want this in units of pixels
    double dl=0.1; //stepsize in uas in x and y
    double dlxy=sqrt(2.)*dl;
    
#pragma omp parallel for schedule(dynamic,1) collapse(2) shared(x,y,vx,vy,imageI,L,dl)
    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){

	    double totl=0.0;
	    //starting point for position and velocity in the pixel center
	    double x_line=x[i];
	    double y_line=y[j];
	    double vx_line=vx[i][j];
	    double vy_line=vy[i][j];
	    double Ilocal;
	    double h;
	    double hsum=0.0;
	    double LP=sqrt(imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2])/imageS[i][j][0];
	    while(totl < 0.5*L && x_line>startx && y_line>starty && x_line<stopx && y_line<stopy){
		x_line += vx_line*dl;
		y_line += vy_line*dl;
		totl += dlxy;
		Ilocal=interp_scalar2d(x_line,y_line,imageI);
		h=1./L;//hanning_windowed(totl,totl+dlxy);
		hsum += h;
		image[i][j] += Ilocal*h;
		vx_line=interp_scalar2d(x_line,y_line,vx);
		vy_line=interp_scalar2d(x_line,y_line,vy);
//		fprintf(stdout,"%g %g \n",x_line,y_line);
	    }
	    totl=0.0;
	    //starting point for position and velocity
	    x_line=x[i];
	    y_line=y[j];
	    vx_line=vx[i][j];
	    vy_line=vy[i][j];
	    while(totl < 0.5*L && x_line>startx && y_line>starty && x_line<stopx && y_line<stopy){
		x_line += -vx_line*dl;
		y_line += -vy_line*dl;
		totl += dlxy;
		Ilocal=interp_scalar2d(x_line,y_line,imageI);
		h=1./L;//hanning_windowed(totl,totl+dlxy);
		hsum+=h;
		image[i][j] += Ilocal*h;
		vx_line=interp_scalar2d(x_line,y_line,vx);
		vy_line=interp_scalar2d(x_line,y_line,vy);
//		fprintf(stdout,"%g %g \n",x_line,y_line);
	    }
	    image[i][j]/=hsum;
//	    if(LP > 0.5){
//		image[i][j]=log(imageS[i][j][0]+1e-50);
//	    }
	    if (j==0) fprintf(stderr,"%d ",i);
	    
	    
	}
    }
    fprintf(stderr,"\n");
    
    fprintf(stderr,"done\n");
    make_ppm(image, freq, "ipole_lic.ppm");

    
    fp = fopen("ipole2.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "unable to open %s\n", fname);
        exit(1);
    }
     for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            fprintf(fp, "%d %d %g %g %g %g %g  %g\n",
                    i, j,
                    image[i][j],
                    image[i][j],
                    imageS[i][j][1],
                    imageS[i][j][2],
                    imageS[i][j][3],
                    1.
		);
	}
	fprintf(fp, "\n");
     }
     fclose(fp);



    
}


//convert x to ij coordinates of a figure
void xtoij(double xi, double yi, int *i, int *j, double *delx, double *dely)
{
    
    //find index of pixel  
    *i = (int) ((xi - startx) / dx - 0.5 + 1000) - 1000;
    *j = (int) ((yi - starty) / dy - 0.5 + 1000) - 1000;

    //and distance from boundaries
    if(*i < 0) {
	*i = 0 ;
	*delx = 0. ;
    } else if(*i > NX-1) {
	*i = NX-1 ;
	*delx = 1. ;
    } else {
	*delx = (xi - ((*i + 0.5) * dx + startx)) / dx;
    }

    
    if(*j < 0) {
	*j = 0 ;
	*dely = 0. ;
    } else if(*j > NY-1) {
	*j = NY-1 ;
	*dely = 1. ;
    } else {
	*dely = (yi - ((*j + 0.5) * dy + starty)) / dy;
    }

    
    return;
    
}

//interpolate vector or scalar to the point
double interp_scalar2d(double xi, double yi, double var[NX][NY])
{
    int i, j, ip1, jp1;
    double delx, dely, b1, b2, interp;

    xtoij(xi, yi, &i, &j, &delx, &dely);
    
    ip1 = i+1;
    jp1 = j+1;
    b1 = 1.-delx;
    b2 = 1.-dely;
    
    interp = var[i][j]*b1*b2 +
	var[i][jp1]*b1*dely +
	var[ip1][j]*delx*b2 +
	var[ip1][jp1]*delx*dely;
        
    return(interp);
    
}

double hanning_windowed(double a,double b){

    double c,d,beta;

    c=1.;// hanning window 
    d=2;//40.;
    beta=0.0;
    
    double h=0.25*(
	b-a+(sin(b*c)-sin(a*c))/c+
	(sin(b*d+beta)-sin(a*d+beta))/d+
	(sin(b*(c-d)-beta)-sin(a*(c-d)-beta))/2./(c-d)+
	(sin(b*(c+d)+beta)-sin(a*(c+d)+beta))/2./(c+d)
	);
	
    return(h);
}
