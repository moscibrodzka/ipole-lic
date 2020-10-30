#include "decs.h"
#include <omp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

double imageS[NX][NY][NDIM];
double image[NX][NY];
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

int main(int argc, char *argv[]){
    
    FILE *fp;
    int i, j;    
    char* fname;
    double dum_dbl;
    int dum_int;
    double freq=230.e9;

    
    fp = fopen(argv[1], "r");
    if (fp == NULL) {
    	fprintf(stderr, "unable to open %s\n", fname);
    	exit(1);
    }

    //ignore 9 lines in the file header 
    char ignore[1024];
    for (i=0;i<10;i++) fgets(ignore, sizeof(ignore), fp);
    
    for (j = 0; j < NY; j++) {
	for (i = 0; i < NX; i++) {
	    //if original ipole image
//   	    fscanf(fp, "%d %d %lf %lf %lf %lf %lf %lf\n",
//	    	    &dum_int, &dum_int, &dum_dbl,&imageS[i][j][0],&imageS[i][j][1],&imageS[i][j][2],&imageS[i][j][3],&dum_dbl);
	    //if fits file dumped to txt
   	    fscanf(fp, "%lf %lf %lf %lf %lf %lf  \n",
		   &dum_dbl, &dum_dbl, &imageS[i][NY-j][0],&imageS[i][NY-j][1],&imageS[i][NY-j][2],&imageS[i][NY-j][3]);
    	}
//    	fscanf(fp,"\n");
    }
    fclose(fp);

    fprintf(stdout,"read %s \n",argv[1]);


    
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
	    image[i][j] = pow(imageS[i][j][0],1.0);
	    //image[i][j] = sqrt(pow(imageS[i][j][1],2.0)+pow(imageS[i][j][2],2.0));
    make_ppm(image, freq, "image_fnu.ppm");
    fprintf(stdout,"made image_fnu.ppm \n");

    
    //in [uas]
    for (int i = 0; i < NX; i++){
	x[i]=-FOV/2. + (i+0.5)*dx;
	y[i]=-FOV/2. + (i+0.5)*dy;	
    }
    
    double amp,ang,Imax,QUmax;
    Imax=0.0;
    QUmax=0.0;
    //construct a vector field for each pixel
    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){
	    amp=1.;//FOV/NX;
	    ang=0.5*atan2(imageS[i][j][2],imageS[i][j][1]);//+M_PI/2.; //here one can att 90 deg to plot BVPAs
	    vx[i][j] = (-sin(ang)*amp);
	    vy[i][j] = ( cos(ang)*amp);
	    image[i][j]=0.0; 
	    //here one has many options for the background intensity
//	    imageI[i][j]=rand();// uniform 
	    imageI[i][j]=imageS[i][j][0]*rand(); // brightness ~ I
//	    imageI[i][j]=imageS[i][j][0]; // brightness ~ I
//	    imageI[i][j]=(imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2])*rand(); // ~P^2
//	    imageI[i][j]=sqrt(imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2])*rand();// ~(P)
//	    imageI[i][j]=rand()*pow(imageS[i][j][0],0.1); //some power of I
	}
    }

    // define length scale on the files L in uas to integrate over
    double L=20.;//4.;//3;//10.;//80.; //uas  //parameter, the larger the number the smoother images
    double dl=0.01; //stepsize in uas in x and y
    double dlxy=sqrt(2.)*dl;
    
#pragma omp parallel for schedule(dynamic,1) collapse(2) shared(x,y,vx,vy,imageI,dl,L)
    for (int i = 0; i < NX-10; i++){
        for (int j = 0; j < NY; j++){
	    double totl1=0.0;
	    //starting point for position and velocity in the pixel center
	    double x_line=x[i];
	    double y_line=y[j];
	    double vx_line=vx[i][j];
	    double vy_line=vy[i][j];
	    double Ilocal;
	    double h;
	    double hsum=0.0;
            //integrate forward until we hit boundary of image or integration limit
	    while(totl1 < 0.5*L){
	    	x_line += vx_line*dl;
	    	y_line += vy_line*dl;
		if(x_line < startx && y_line < starty && x_line > stopx && y_line > stopy) break;
	    	totl1 += dlxy;
		Ilocal=interp_scalar2d(x_line,y_line,imageI);
		h=1./L;
	    	//here one can use another filter
	    	//h=hanning_windowed(totl,totl+dlxy);
		hsum += h;
		image[i][j] += Ilocal*h;
		vx_line=interp_scalar2d(x_line,y_line,vx);
	    	vy_line=interp_scalar2d(x_line,y_line,vy);
	    }
	    double totl2=0.0;
	    //starting point for position and velocity
	    x_line=x[i];
	    y_line=y[j];
	    vx_line=vx[i][j];
	    vy_line=vy[i][j];
	    while(totl2 < 0.5*L){
	    	x_line += -vx_line*dl;
	    	y_line += -vy_line*dl;
		if(x_line < startx && y_line < starty && x_line > stopx && y_line > stopy) break;
	    	totl2 += dlxy;
	    	Ilocal=interp_scalar2d(x_line,y_line,imageI);
	    	h=1./L;
	    	//here one can use another filter
	    	//h=hanning_windowed(totl,totl+dlxy);
		hsum+=h;
		image[i][j] += Ilocal*h;
	    	vx_line=interp_scalar2d(x_line,y_line,vx);
	    	vy_line=interp_scalar2d(x_line,y_line,vy);
	    }
	    image[i][j]/=hsum;
	}
    }

    //normalize both images to add them up and somewhat conserve the flux
    double image_max=0.0;
    double image_max2=0.0;
    double QU_max2=0.0;
    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){
	    if(image[i][j] > image_max) image_max=image[i][j];
	    if(imageS[i][j][0] > image_max2) image_max2=imageS[i][j][0];
	    if(imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2] > QU_max2){
		QU_max2=imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2];
	    }
	}
    }

    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){
	    image[i][j]=image[i][j]/image_max;
	    imageS[i][j][0]=imageS[i][j][0]/image_max2;
	}
    }

    //add weight of stroke
    double fac;
    for (int i = 0; i < NX; i++){
        for (int j = 0; j < NY; j++){
	    //press brush hard only in regions with high pol flux
	    fac=(imageS[i][j][1]*imageS[i][j][1]+imageS[i][j][2]*imageS[i][j][2])/(QU_max2);
	    image[i][j] = imageS[i][j][0]*(1.-fac) + image[i][j]*fac;
	}
    }
    
    
    make_ppm(image, freq, "image_lic.ppm");
    fprintf(stdout,"made image_lic.ppm \n");
    
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
