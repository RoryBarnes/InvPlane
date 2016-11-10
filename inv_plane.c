/**********************************************************************
 *                                                                    *
 * INV_PLANE.C                                                        *
 *                                                                    *
 * Author: Rory Barnes (rory@astro.washington.edu)                    *
 *                                                                    *
 * To compile: gcc -o invplane inv_plane.c -lm                        *
 *                                                                    *
 * This code tranforms a system in which inclinations and longitudes  *
 * of ascending node is measured from an arbitrary reference frame    *
 * into one in which they are measured from the invariable plane.     *
 * The input file  contains N+2 lines, where N is the number of       *
 * orbiters (don't count the primary!). The first line is N, the      *
 * second line is the mass of the central body, and the next N lines  *
 * contain the orbital parameters of the orbiters in the format:      *
 * Mass SemimajorAxis Eccentricity Inclination LongitudeAscendingNode *
 * ArgumentPericenter MeanAnomaly. The units are solar masses, AU,    *
 * and degrees.                                                       *
 *                                                                    *
 * The user specifies the output format at the command line:          *
 * -a   an ASCII text file                                            *
 * -m   a big.in file for use with MERCURY                            *
 * -h   a .hnb file for use with HNBODY                               *
 * -v   report all warnings                                           *
 *                                                                    *
 * For example, to transform a system with parameters provided in a   *
 * file called system.in into a big.in file for use with MERCURY,     *
 * the command is: invplane -m system.in.                             *
 *                                                                    *
 * Thanks to Russell Deitrick and Pramod Gupta for testing.           *
 *                                                                    *
 **********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define dot(a,b)        (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define cross_1(z,a,b)  z[0]= a[1]*b[2]-a[2]*b[1]
#define cross_2(z,a,b)  z[1]= a[2]*b[0]-a[0]*b[2]
#define cross_3(z,a,b)  z[2]= a[0]*b[1]-a[1]*b[0]
#define cross(z,a,b)    cross_1(z,a,b);cross_2(z,a,b);cross_3(z,a,b)
#define MJUP  1.8987e30
#define MEARTH 5.9742e27
#define MSUN 1.98892e33
#define AUCM 1.49598e13
#define BIGG 6.672e-8
#define EPS 1e-2       // Limit to produce warnings

typedef struct elem_struct
{
  double a;       // Semi-major axis
  double e;       // Eccentricity
  double i;       // Inclination
  double lasc;    // Longitude of ascending node
  double aper;    // Arg. of Pericenter
  double mean_an; // Mean Anomaly
} ELEMS;

/* Compute orbital elements from Cartesian coordinates */
void elems(double mu, double *x,double *v,ELEMS *elem) {
   double hx,hy,hz,hsq,hxy,h,r,vsq,rdot;
   double sin_lasc,cos_lasc,sin_aperf,cos_aperf,sinf,cosf;
   double a,e,sin_ecc,cos_ecc,arg;
   int i;
/* 
 * Compute various intermediate quantities.
 */
   hx= (double)(x[1]*v[2]-x[2]*v[1]);
   hy= (double)(x[2]*v[0]-x[0]*v[2]);
   hz= (double)(x[0]*v[1]-x[1]*v[0]);
   hsq= hx*hx+hy*hy+hz*hz;
   h= sqrt(hsq);
   hxy= sqrt(hx*hx+hy*hy);
   r= sqrt((double)dot(x,x));
   vsq= (double)dot(v,v);
   rdot= (double)dot(x,v);
   a= 1/(2/r-vsq/mu);
   e= sqrt(1-hsq/(mu*a));

   sin_lasc= hx/hxy;
   cos_lasc= -hy/hxy;
   sin_aperf= (double)x[2]*h/(hxy*r);
   cos_aperf= ((double)x[0]*cos_lasc+(double)x[1]*sin_lasc)/r;
   cosf= (hsq/(mu*r)-1)/e;
   if(fabs(cosf)>1) {
     if(cosf>0) 
       cosf= 1;
     else
       cosf= -1; 
   }
   sinf= sqrt(1-cosf*cosf);
   if(rdot<0)  
     sinf= -sinf;
   cos_ecc= (1-r/a)/e;
   if (fabs(cos_ecc)>1) {
      if(cos_ecc>0)
        cos_ecc= 1;
      else
        cos_ecc= -1;
   }
   sin_ecc = sqrt(1-cos_ecc*cos_ecc);
   if (rdot<0) 
     sin_ecc = -sin_ecc;

   elem->a=a;
   elem->e=e;
   elem->i=atan2(hxy/h,hz/h);
   elem->lasc= atan2(sin_lasc,cos_lasc);
   arg= atan2(sin_aperf,cos_aperf)-atan2(sinf,cosf);
   while(arg>M_PI) 
     arg-= 2*M_PI;
   while(arg<-M_PI) 
     arg+= 2*M_PI;
   elem->aper = arg;
   elem->mean_an = atan2(sin_ecc,cos_ecc)-e*sin_ecc;
}

double distance(double *x)  {
  return sqrt(dot(x,x));
}

/* Compute Cartesian coordinates from orbital elements */
void cartes(double *x,double *v,ELEMS *elem,double mu) {
  double a,e,m,cosi,sini,cos_lasc,sin_lasc,cos_aper,sin_aper;
  double es,ec,w,wp,wpp,wppp,ecc,dx,lo,up,next;
  int iter;
  double sin_ecc,cos_ecc,l1,m1,n1,l2,m2,n2;
  double xi,eta,vel_scl;
 
  a=elem->a;
  e=elem->e;  
  m=elem->mean_an;
  cosi=cos(elem->i);
  sini=sin(elem->i);
  cos_lasc=cos(elem->lasc);
  sin_lasc=sin(elem->lasc);
  cos_aper=cos(elem->aper);
  sin_aper=sin(elem->aper);
 
  /*
   * Reduce mean anomoly to [0, 2*PI) 
   */
  m-= ((int)(m/(2*M_PI)))*2*M_PI;
  /*
   * Solve Kepler's equation.
   */
  if(sin(m)>0)
    ecc= m+0.85*e;
  else 
    ecc= m-0.85*e;
  lo=-2*M_PI;
  up= 2*M_PI;

  for(iter= 1;iter<=32;iter++) {   
    es=e*sin(ecc);
    ec=e*cos(ecc);
    w=ecc-es-m;
    wp=1-ec;
    wpp=es;
    wppp=ec;
    if(w>0)
      up=ecc;
    else 
      lo= ecc;
    dx=-w/wp;
    dx=-w/(wp+dx*wpp/2);
    dx=-w/(wp+dx*wpp/2+dx*dx*wppp/6);
    next=ecc+dx;
    if(ecc==next)
        break;
    if((next>lo) && (next<up))
      ecc=next;
    else 
      ecc= (lo+up)/2;
    if((ecc==lo)||(ecc==up))
        break;
    if(iter>30)
      printf("%4d %23.20f %e\n",iter,ecc,up-lo);
  }

  if(iter>32) {
    printf("Kepler soln failed.\n");
    exit(1);
    }
    
  cos_ecc=cos(ecc);
  sin_ecc=sin(ecc);
    
  l1=cos_lasc*cos_aper-sin_lasc*sin_aper*cosi;
  m1=sin_lasc*cos_aper+cos_lasc*sin_aper*cosi;
  n1=sin_aper*sini;
  l2=-cos_lasc*sin_aper-sin_lasc*cos_aper*cosi;
  m2=-sin_lasc*sin_aper+cos_lasc*cos_aper*cosi;
  n2=cos_aper*sini;
    
  xi=a*(cos_ecc-e);eta= a*pow(1-e*e,0.5)*sin_ecc;
  x[0]=l1*xi+l2*eta;
  x[1]=m1*xi+m2*eta;
  x[2]=n1*xi+n2*eta;
  vel_scl=pow((mu*a)/dot(x,x),0.5);
  xi=-vel_scl*sin_ecc;eta= vel_scl*pow(1-e*e,0.5)*cos_ecc;
  v[0]=l1*xi+l2*eta;
  v[1]=m1*xi+m2*eta;
  v[2]=n1*xi+n2*eta;
}

/* Astrocentric -> Barycentric Cartesian coordinates */
void hel_bar(double **hel,double **bar,double *m,double *ms,int P)
{
   int i,p;

   for(i= 0;i<3;i++)
      bar[0][i] = 0;
   for(p= 1;p<=P;p++)
      for(i=0;i<3;i++)
         bar[0][i] -= m[p]/ms[P]*hel[p][i];
   for(p= 1;p<=P;p++)
      for(i=0;i<3;i++)
         bar[p][i] = hel[p][i] + bar[0][i];
}

/* Barycentric -> Astrocentric Cartesian coordinates */
void bar_hel(double **bar,double **hel,int P)
{
   int i,p;

   for (p=1;p<=P;p++)
      for (i=0;i<3;i++) hel[p][i] = bar[p][i] - bar[0][i];
   for (i=0;i<3;i++) hel[0][i] = 0;

}

/* Calculate the total angular momentum vector */
double* invariable_plane(double **x,double **v,int P,double *m) {
  int p,i;
  double *dh,*h;
  double mg;

  dh=malloc(3*sizeof(double));
  h=malloc(3*sizeof(double));

  cross(h,x[0],v[0]);
  for (i=0;i<3;i++)
    h[i] *= m[0];   
  for(p=1;p<=P;p++) {
    cross(dh,x[p],v[p]);
    for(i=0;i<3;i++)
      h[i]+= m[p]*dh[i];
  }
  mg= pow(dot(h,h),0.5);
  for(i=0;i<3;i++) 
    h[i]/= mg;
  return h;
}

/* Rotate coordinates */
void rotate(double *z,double **x,int np)  {
    double phi, theta;
    int i,k;
    double **x1;
    
    x1=malloc((np+1)*sizeof(double*));
    for (i=0;i<=np;i++)
        x1[i]=malloc(3*sizeof(double));

    theta = atan2(z[1],z[0]);
    phi = atan2(sqrt(z[0]*z[0] + z[1]*z[1]),z[2]);

    /* Rotate about z-axis */
    for (i=0;i<=np;i++) {
      x1[i][0] = x[i][0]*cos(theta) + x[i][1]*sin(theta);
      x1[i][1] = -x[i][0]*sin(theta) + x[i][1]*cos(theta);
      x1[i][2] = x[i][2];
    }
    
    /* Rotate about new y-axis (z -> x, x -> y) */
    for (i=0;i<=np;i++) {
      x[i][0] = -x1[i][2]*sin(phi) + x1[i][0]*cos(phi);
      x[i][1] = x1[i][1];
      x[i][2] = x1[i][2]*cos(phi) + x1[i][0]*sin(phi);
    }

    free(x1);
}

int main(int argc, char *argv[]) {
  int k,j,np;
  ELEMS *p;
  double *m,*ms,*mu;
  double **hex,**hev,**bax,**bav;
  double *zprime,*d0,*df,*vmag0,*vmagf,ksq;
  double *dist0,*distf,*phi0,*phif;
  int c,dobig=0,dohnb=0,doascii=0,verbose=0;  

  char id[16],outfile[256];
  FILE *ifp,*afp,*mfp,*hfp;
  
  if (argc > 6) {
      (void) fprintf(stderr,"Usage: %s -mhav file\n",argv[0]);
      exit(1);
  }
  /* This could be done better. The user could supply unknown flag and 
     nothing happens! */
  for (k=0;k<argc;k++) {
    if (strcmp(argv[k],"-m") == 0) 
      dobig=1;
    if (strcmp(argv[k],"-h") == 0)
      dohnb=1;
    if (strcmp(argv[k],"-a") == 0)
      doascii=1;
    if (strcmp(argv[k],"-v") == 0)
      verbose=1;
  }

  if (dobig==0 && dohnb==0 && doascii==0) {
    fprintf(stderr,"ERROR: Must specify output file type: Mercury (-m), HNBody (-h) or ascii (-a)\n");
    exit(1);
  }

  ifp=fopen(argv[argc-1],"r");

  fscanf(ifp,"%d",&np);
  /* Initialize arrays */
  m=malloc((np+1)*sizeof(double));
  ms=malloc((np+1)*sizeof(double));
  mu=malloc((np+1)*sizeof(double));
  hex=malloc((np+1)*sizeof(double*));
  hev=malloc((np+1)*sizeof(double*));
  bax=malloc((np+1)*sizeof(double*));
  bav=malloc((np+1)*sizeof(double*));
  p=malloc((np+1)*sizeof(ELEMS));
  vmag0=malloc(np*sizeof(double));
  vmagf=malloc(np*sizeof(double));
  for (k=0;k<=np;k++) {
    hex[k]=malloc(3*sizeof(double));
    hev[k]=malloc(3*sizeof(double));
    bax[k]=malloc(3*sizeof(double));
    bav[k]=malloc(3*sizeof(double));

    /* Initialize values. */
    m[j]=0;
    p[j].a=0;
    p[j].e=0;
    p[j].i=0;
    p[j].lasc=0;
    p[j].aper=0;
    p[j].mean_an=0;
  }
  zprime=malloc(3*sizeof(double));
  d0=malloc(np*sizeof(double));
  df=malloc(np*sizeof(double));
  dist0=malloc(np*sizeof(double));
  distf=malloc(np*sizeof(double));
  phi0=malloc(np*sizeof(double));
  phif=malloc(np*sizeof(double));

  fscanf(ifp,"%lf",&m[0]);
  m[0] *= MSUN;

  /* Should add control to prevent a user from inserting more or less than np
     planets. */
      
  for (j=1;j<=np;j++) {   
    /* Read in astrocentric orbital elements */
    /* Mass in Solar units */
    fscanf(ifp,"%lf %lf %lf %lf %lf %lf %lf",&m[j],&p[j].a,&p[j].e,&p[j].i,&p[j].lasc,&p[j].aper,&p[j].mean_an);
    p[j].i *= M_PI/180;
    p[j].lasc *= M_PI/180;
    p[j].mean_an *= M_PI/180;
    p[j].aper *= M_PI/180;
    p[j].a *= AUCM;
    m[j] *= MSUN;
    if (m[j] <= 0) {
      fprintf(stderr,"ERROR: Mass of planet %d must be positive.\n",j);
      exit(1);
    }
    if (p[j].a <= 0) {
      fprintf(stderr,"ERROR: Semi-major axis of planet %d must be positive.\n",j);
      exit(1);
    }
    if (p[j].e < 0) {
      fprintf(stderr,"ERROR: Eccentricity of planet %d cannot be negative.\n",j);
      exit(1);
    }
    if (p[j].e >= 1) {
      fprintf(stderr,"ERROR: Eccentricity of planet %d must be less than 1.\n",j);
      exit(1);
    }
    if (p[j].i >= M_PI) {
      fprintf(stderr,"ERROR: Inclination of planet %d must be less than 180.\n",j);
      exit(1);
    }
    if (p[j].i < 0) {
      fprintf(stderr,"ERROR: Inclination of planet %d must be positive.\n",j);
      exit(1);
    }
  }

  ksq=BIGG*m[0];

  ms[0]=m[0];
  for (k=1;k<=np;k++) {
      ms[k]=ms[k-1]+m[k];
      mu[k]=(ksq*ms[k])/ms[k-1];
  }     

  for (j=1;j<np;j++) {
      phi0[j-1]=acos(cos(p[j].i)*cos(p[j+1].i) + sin(p[j].i)*sin(p[j+1].i)*cos(p[j].lasc - p[j+1].lasc));
      }

  /* Convert to Cartesian astrocentric coordinates */
  for (j=1;j<=np;j++) {
    cartes(hex[j],hev[j],&p[j],mu[1]);
    /* Record initial speed to check final coordinates later */
    vmag0[j-1] = dot(hev[j],hev[j]);
  }

  /* Calculate relative distances between bodies to check later */
  for (j=1;j<np;j++) {
    dist0[j-1]=sqrt(pow((hex[j+1][0]-hex[j][0]),2) + pow((hex[j+1][1]-hex[j][1]),2) + pow((hex[j+1][2]-hex[j][2]),2));
  }

  /* Get distance from origin*/
  for (j=1;j<=np;j++) 
    d0[j-1] = distance(hex[j]);

  /* Initialize stellar position and velocity */
  for (j=0;j<3;j++) {
    hex[0][j]=0.0;
    hev[0][j]=0.0;
  }
    
  /* Convert to Cartesian barycentric coordinates */
  hel_bar(hex,bax,m,ms,np);
  hel_bar(hev,bav,m,ms,np);
  
  /* Calculate invariable plane */
  zprime = invariable_plane(bax,bav,np,m);

  /* Rotate coordinates */
  rotate(zprime,bax,np);
  rotate(zprime,bav,np);
  
  /* Confirm that L is || to z-hat */
  zprime = invariable_plane(bax,bav,np,m);
  if (zprime[2] != 1 || verbose) 
    printf("WARNING: zprime[2] = %lf != 1\n",zprime[2]);

  /* Convert to Cartesian astrocentric coordinates */
  bar_hel(bax,hex,np);
  bar_hel(bav,hev,np);

  /* Get distances, velocities, and assert */
  for (j=1;j<=np;j++) {
    df[j-1]=distance(hex[j]);
    if (fabs((d0[j-1] - df[j-1])/d0[j-1]) > EPS || verbose) 
      printf("WARNING: d0[%d] = %lf, but df[%d] = %lf\n",j-1,d0[j-1],j-1,df[j-1]);
  }

  for (j=1;j<np;j++) {
    distf[j-1]=sqrt(pow((hex[j+1][0]-hex[j][0]),2) + pow((hex[j+1][1]-hex[j][1]),2) + pow((hex[j+1][2]-hex[j][2]),2));
    if (fabs((dist0[j-1] - distf[j-1])/dist0[j-1]) > EPS || verbose) 
      printf("WARNING: dist0[%d] = %lf, but distf[%d] = %lf\n",j-1,dist0[j-1],j-1,distf[j-1]);
    vmagf[j-1]=dot(hev[j],hev[j]);
    if (fabs((vmag0[j-1] - vmagf[j-1])/vmag0[j-1]) > EPS || verbose) 
      printf("WARNING: vmag0[%d] = %lf, but vmagf[%d] = %lf\n",j-1,vmag0[j-1],j-1,vmagf[j-1]);
    }
    
  /* Convert to astrocentric elements */
  for (j=1;j<=np;j++)
    elems(mu[1],hex[j],hev[j],&p[j]);
   
  /* Assert that mutual inclinations are the same */
  for (j=1;j<np;j++) {
    phif[j-1]=acos(cos(p[j].i)*cos(p[j+1].i) + sin(p[j].i)*sin(p[j+1].i)*cos(p[j].lasc - p[j+1].lasc));
    if (fabs((phi0[j-1] - phif[j-1])/phi0[j-1]) > EPS || verbose) 
      printf("WARNING: phi0[%d] = %lf, but phif[%d] = %lf\n",j-1,phi0[j-1],j-1,phif[j-1]);
  }

  /* Convert elements to useful units */
  for (j=1;j<=np;j++) {
    p[j].i *= 180/M_PI;
    while (p[j].i < 0) 
      p[j].i += 360;
    while (p[j].i >= 360) 
      p[j].i -= 360;
    p[j].aper *= 180/M_PI;
    while (p[j].aper < 0) 
      p[j].aper += 360;
    while (p[j].aper >= 360) 
      p[j].aper -= 360;
    p[j].lasc *= 180/M_PI;
    while (p[j].lasc < 0) 
      p[j].lasc += 360;
    while (p[j].lasc >= 360) 
      p[j].lasc -= 360;
    p[j].mean_an *= 180/M_PI;
    while (p[j].mean_an < 0) 
      p[j].mean_an += 360;
    while (p[j].mean_an >= 360) 
      p[j].mean_an -= 360;
    p[j].a /= AUCM;
    m[j] /= MSUN;
  }    

  /* Write ASCII text file in same format as input file */
  if (doascii) {
    sprintf(outfile,"%s.inv",argv[argc-1]);
    afp=fopen(outfile,"w");
    fprintf(afp,"%d\n",np);
    fprintf(afp,"%lf\n",m[0]/MSUN);
    for (j=1;j<=np;j++) 
      fprintf(afp,"%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n",m[j],p[j].a,p[j].e,p[j].i,p[j].lasc,p[j].aper,p[j].mean_an);
    if (verbose) 
      printf("Wrote %s.\n",outfile);
  }

  /* Write big.in file for use with MERCURY */
  if (dobig) {
    mfp=fopen("big.in","w");
    fprintf(mfp,")O+_06 Rotated coordinates for file %s\n",argv[argc-1]);
    fprintf(mfp,")---------------------------------------------------------------------\n");
    fprintf(mfp, " style (Cartesian, Asteroidal, Cometary) = Asteroidal\n");
    fprintf(mfp,"  epoch (in days) = 0.0\n");
    fprintf(mfp,")---------------------------------------------------------------------\n");
    for (j=1;j<=np;j++) {
      fprintf(mfp," %d m=%.5e\n",j,m[j]);
      fprintf(mfp," %.5e %.5e %.5e\n",p[j].a,p[j].e,p[j].i);
      fprintf(mfp," %.5e %.5e %.5e\n",p[j].aper,p[j].lasc,p[j].mean_an);
      fprintf(mfp," 0. 0. 0.\n");
    }
    if (verbose) 
      printf("Wrote big.in.\n");
  }

  /* Write HNB file for use with HNBody */
  if (dohnb) {
    sprintf(outfile,"%s.hnb",argv[argc-1]);
    hfp=fopen(outfile,"w");
    fprintf(hfp,"AngleUnit:  deg\n");
    fprintf(hfp,"LengthUnit: AU\n");
    fprintf(hfp,"MassUnit:   Msun\n");
    fprintf(hfp,"TimeUnit:   yr\n");
    fprintf(hfp,"StepSize:   \n");
    fprintf(hfp,"M  =  %lf\n",m[0]/MSUN);
    fprintf(hfp,"N  = %d\n\n",np+1);
    fprintf(hfp,"InputOrder: Mass SemiMajorAxis  Eccentricity Inclination \\\n");
    fprintf(hfp,"\tLongAscendNode ArgPeriapse  MeanAnomaly\n");
    fprintf(hfp,"ParticleType: HWP\n");
    for (j=1;j<=np;j++) 
      fprintf(hfp,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",m[j],p[j].a,p[j].e,p[j].i,p[j].lasc,p[j].aper,p[j].mean_an);
    
    fprintf(hfp,"\nTfinal: \n");
    fprintf(hfp,"OutputFiles: %%d.dat\n");
    fprintf(hfp,"OutputOrder: time Semi Ecc Inc LongA ArgP MeanAn\n");          
    fprintf(hfp,"OutputInterval: \n");
    fprintf(hfp,"OutputCoord: Bodycentric\n");
    fprintf(hfp,"OutputHeader: True\n");
    fprintf(hfp,"OutputDigits: 8\n");
    fprintf(hfp,"OutputTypes: HWPs\n");
    if (verbose) 
      printf("Wrote %s.\n",outfile);
  }

    return 0;
}
