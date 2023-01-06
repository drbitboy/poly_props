
/* Update 2023-01-06 by Brian T. Carcich, Latchmoor Services, INC.
 *
 * Build:  gcc volInt.c -o volInt -lm
 *
 * Updated URL for Brian Mirtich:
 *
 *   https://people.eecs.berkeley.edu/~jfc/mirtich/
 *
 * Updated URL for this app:
 *
 *   https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html
 *
 */

        /*******************************************************
        *                                                      *
        *  volInt.c                                            *
        *                                                      *
        *  This code computes volume integrals needed for      *
        *  determining mass properties of polyhedral bodies.   *
        *                                                      *
        *  For more information, see the accompanying README   *
        *  file, and the paper                                 *
        *                                                      *
        *  Brian Mirtich, "Fast and Accurate Computation of    *
        *  Polyhedral Mass Properties," journal of graphics    *
        *  tools, volume 1, number 1, 1996.                    *
        *                                                      *
        *  This source code is public domain, and may be used  *
        *  in any way, shape or form, free of charge.          *
        *                                                      *
        *  Copyright 1995 by Brian Mirtich                     *
        *                                                      *
        *  mirtich@cs.berkeley.edu                             *
        *  http://www.cs.berkeley.edu/~mirtich                 *
        *                                                      *
        *******************************************************/

/*
        Revision history

        26 Jan 1996     Program creation.

         3 Aug 1996     Corrected bug arising when polyhedron density
                        is not 1.0.  Changes confined to function main().
                        Thanks to Zoran Popovic for catching this one.

        27 May 1997     Corrected sign error in translation of inertia
                        product terms to center of mass frame.  Changes
                        confined to function main().  Thanks to
                        Chris Hecker.

        06 Jan 2023     Updated URLs, replace tabs with spaces, remove
                        trailing whitespace, move to Github,
                        by B.T. Carcich
*/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

/*
   ============================================================================
   constants
   ============================================================================
*/
#ifdef USEMAX_
#undef USEMAX_
#endif

//#define USEMAX_           // uncomment to used fixed-size arrays

#ifdef USEMAX_
#  define MAX_VERTS 100     /* maximum number of polyhedral vertices */
#  define MAX_FACES 100     /* maximum number of polyhedral faces */
#endif

#define MX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */

#define X 0
#define Y 1
#define Z 2

/*
   ============================================================================
   macros
   ============================================================================
*/

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/*
   ============================================================================
   data structures
   ============================================================================
*/

typedef struct {
  int numVerts;
  double norm[3];
  double w;
  int iverts[MX_POLYGON_SZ];
  struct polyhedron *poly;
} FACE;

typedef struct polyhedron {
  int numVerts, numFaces;
# ifdef USEMAX_
    double verts[MAX_VERTS][3];
    FACE faces[MAX_FACES];
# else
    double **verts;
    FACE *faces;
# endif
} POLYHEDRON;


/*
   ============================================================================
   globals
   ============================================================================
*/

static int A;   /* alpha */
static int B;   /* beta */
static int C;   /* gamma */

/* projection integrals */
static double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
static double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
static double T0, T1[3], T2[3], TP[3];

double VolIntBinSearchCubicRoot(double CubicCoeffs[4], double lox, double hix)
{
#define CC3 CubicCoeffs[3]
#define CC2 CubicCoeffs[2]
#define CC1 CubicCoeffs[1]
#define CC0 CubicCoeffs[0]
#define EVALCUBIC(X) (CC0+(X)*(CC1+(X)*(CC2+(X)*CC3)))
#define PRINTFCC printf( "CubicCoeffs[3,2,1,0]=\n  %lg\n  %lg\n  %lg\n  %lg\n", CC3, CC2, CC1, CC0)

double lo, hi, mid;
double midx = (lox + hix) / 2.0;
long iter = 0;

  if ( (lo = EVALCUBIC(lox)) == 0.0) return lox;
  if ( (hi = EVALCUBIC(hix)) == 0.0) return hix;

  if ( lo > 0.0 || hi < 0.0)
  {
    printf( "lox,lo=%35.25lg,%35.25lg\nhix,hi=%35.25lg,%35.25lg\n", lox, lo, hix, hi);
    PRINTFCC;
    perror( "VolIntBinSearchCubicRoot:  Bad inputs\n");
    return  midx;
  }

  while ( lox != midx && midx != hix)
  {
    if ( ++iter > 1000)
    {
      printf( "lox,lo=%35.25lg,%35.25lg\nhix,hi=%35.25lg,%35.25lg\n", lox, lo, hix, hi);
      PRINTFCC;
      perror( "VolIntBinSearchCubicRoot:  exceeded 1000 iterations\n");
      return  midx;
    }
    mid = EVALCUBIC(midx);
    if (mid == 0.0) break;
    if (mid > 0.0) hix = midx; else lox = midx;
    midx = (lox + hix) / 2.0;

#   ifdef _DEBUG
    printf( ".");
#   endif

  }

# ifdef _DEBUG
  if ( iter) printf( "\n");
# endif

  return midx;
}

void VolIntCubicRoots(double CubicCoeffs[4], double Roots[3])
{
double inflectx = CC2 / 3.0;
double determx = SQR(inflectx) + (CC1 / 3.0);
double loxmincub;
double hixmaxcub;

  if (determx < 0.0)
  {
    PRINTFCC;
    perror( "VolIntCubicRoots:  negative determinant\n");
    return;
  }
  determx = sqrt(determx);

  if ( determx == 0.0)
  {
    Roots[0] = Roots[1] = Roots[2] = inflectx;
    return;
  }
  loxmincub = inflectx - determx;
  hixmaxcub = inflectx + determx;

# ifdef _DEBUG
  printf( "inflectx  = %35.25lg\n", inflectx);
  printf( "determx   = %35.25lg\n", determx);
  printf( "loxmincub = %35.25lg\n", loxmincub);
  printf( "hixmaxcub = %35.25lg\n", hixmaxcub);
  printf( "cubcoeffs = %lg %lg %lg %lg\n", CC3, CC2, CC1, CC0);
  printf( "\n");
# endif

  Roots[0] = VolIntBinSearchCubicRoot(CubicCoeffs, loxmincub, loxmincub-2.0*determx);
  Roots[1] = VolIntBinSearchCubicRoot(CubicCoeffs, loxmincub, hixmaxcub);
  Roots[2] = VolIntBinSearchCubicRoot(CubicCoeffs, hixmaxcub+2*determx, hixmaxcub);

  return;
}

/*
   ============================================================================
   read in a polyhedron
   ============================================================================
*/

void readPolyhedron(char *name, POLYHEDRON *p)
{

  FILE *fp;
  //char line[200], *c;
  int i, j;  //, n;
  double dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len;
  FACE *f;

  if ( !name) {
    fp = stdin;
  } else if (!(fp = fopen(name, "r"))) {
    printf("i/o error\n");
    exit(1);
  }

  fscanf(fp, "%d", &p->numVerts);
  printf("Reading in %d vertices\n", p->numVerts);


# ifndef USEMAX_
  {
  double *pVerts = (double *) malloc( p->numVerts * ((3 * sizeof(double)) + sizeof(double *)) );
  double **ppVerts = (double **) (pVerts +  (p->numVerts * 3));
    p->verts = ppVerts;
    while ( ppVerts < (p->verts + p->numVerts) )
    {
      *(ppVerts++) = pVerts;
      pVerts += 3;
    }
  }
# endif

  for (i = 0; i < p->numVerts; i++)
    fscanf(fp, "%lf %lf %lf"
          , &p->verts[i][X], &p->verts[i][Y], &p->verts[i][Z]);

  fscanf(fp, "%d", &p->numFaces);
  printf("Reading in %d faces\n", p->numFaces);

# ifndef USEMAX_
  p->faces = (FACE *) malloc( p->numFaces * sizeof(FACE) );
# endif

  for (i = 0; i < p->numFaces; i++) {
    f = &p->faces[i];
    f->poly = p;
    fscanf(fp, "%d", &f->numVerts);
    for (j = 0; j < f->numVerts; j++) fscanf(fp, "%d", &f->iverts[j]);

    /* compute face normal and offset w from first 3 vertices */
    dx1 = p->verts[f->iverts[1]][X] - p->verts[f->iverts[0]][X];
    dy1 = p->verts[f->iverts[1]][Y] - p->verts[f->iverts[0]][Y];
    dz1 = p->verts[f->iverts[1]][Z] - p->verts[f->iverts[0]][Z];
    dx2 = p->verts[f->iverts[2]][X] - p->verts[f->iverts[1]][X];
    dy2 = p->verts[f->iverts[2]][Y] - p->verts[f->iverts[1]][Y];
    dz2 = p->verts[f->iverts[2]][Z] - p->verts[f->iverts[1]][Z];
    nx = dy1 * dz2 - dy2 * dz1;
    ny = dz1 * dx2 - dz2 * dx1;
    nz = dx1 * dy2 - dx2 * dy1;
    len = sqrt(nx * nx + ny * ny + nz * nz);
    f->norm[X] = nx / len;
    f->norm[Y] = ny / len;
    f->norm[Z] = nz / len;
    f->w = - f->norm[X] * p->verts[f->iverts[0]][X]
           - f->norm[Y] * p->verts[f->iverts[0]][Y]
           - f->norm[Z] * p->verts[f->iverts[0]][Z];

  }

  if (name) fclose(fp);

}

/*
   ============================================================================
   compute mass properties
   ============================================================================
*/


/* compute various integrations over projection of face */
void compProjectionIntegrals(FACE *f)
{
  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int i;

  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

  for (i = 0; i < f->numVerts; i++) {
    a0 = f->poly->verts[f->iverts[i]][A];
    b0 = f->poly->verts[f->iverts[i]][B];
    a1 = f->poly->verts[f->iverts[(i+1) % f->numVerts]][A];
    b1 = f->poly->verts[f->iverts[(i+1) % f->numVerts]][B];
    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
    a1_2 = a1 * a1; a1_3 = a1_2 * a1;
    b1_2 = b1 * b1; b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

    P1 += db*C1;
    Pa += db*Ca;
    Paa += db*Caa;
    Paaa += db*Caaa;
    Pb += da*Cb;
    Pbb += da*Cbb;
    Pbbb += da*Cbbb;
    Pab += db*(b1*Cab + b0*Kab);
    Paab += db*(b1*Caab + b0*Kaab);
    Pabb += da*(a1*Cabb + a0*Kabb);
  }

  P1 /= 2.0;
  Pa /= 6.0;
  Paa /= 12.0;
  Paaa /= 20.0;
  Pb /= -6.0;
  Pbb /= -12.0;
  Pbbb /= -20.0;
  Pab /= 24.0;
  Paab /= 60.0;
  Pabb /= -60.0;
}

void
compFaceIntegrals(FACE *f)
{
  double *n, w;
  double k1, k2, k3, k4;

  compProjectionIntegrals(f);

  w = f->w;
  n = f->norm;
  k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

  Fa = k1 * Pa;
  Fb = k1 * Pb;
  Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

  Faa = k1 * Paa;
  Fbb = k1 * Pbb;
  Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
        + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faaa = k1 * Paaa;
  Fbbb = k1 * Pbbb;
  Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab
         + 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
         + 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
         + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faab = k1 * Paab;
  Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
  Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
         + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
}

void compVolumeIntegrals(POLYHEDRON *p)
{
  FACE *f;
  double nx, ny, nz;
  int i;

  T0 = T1[X] = T1[Y] = T1[Z]
     = T2[X] = T2[Y] = T2[Z]
     = TP[X] = TP[Y] = TP[Z] = 0;

  for (i = 0; i < p->numFaces; i++) {

    f = &p->faces[i];

    nx = fabs(f->norm[X]);
    ny = fabs(f->norm[Y]);
    nz = fabs(f->norm[Z]);
    if (nx > ny && nx > nz) C = X;
    else C = (ny > nz) ? Y : Z;
    A = (C + 1) % 3;
    B = (A + 1) % 3;

    compFaceIntegrals(f);

    T0 += f->norm[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));

    T1[A] += f->norm[A] * Faa;
    T1[B] += f->norm[B] * Fbb;
    T1[C] += f->norm[C] * Fcc;
    T2[A] += f->norm[A] * Faaa;
    T2[B] += f->norm[B] * Fbbb;
    T2[C] += f->norm[C] * Fccc;
    TP[A] += f->norm[A] * Faab;
    TP[B] += f->norm[B] * Fbbc;
    TP[C] += f->norm[C] * Fcca;
  }

  T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
  T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
  TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;
}


/*
   ============================================================================
   main
   ============================================================================
*/


int main(int argc, char *argv[])
{
double eigVals[3];      /* eigenvalues */
double eigVecs[3][3];   /* eigen vectors */
int iEig;
char fn[] = { "cube.txt" };
#ifdef _DEBUG
#define FN fn
#else
#define FN 0
#endif

  POLYHEDRON p;
  double density, mass;
  double r[3];            /* center of mass */
  double J[3][3];         /* inertia tensor */

  double CubicCoeffs[4];  /* cubic coefficients */

  readPolyhedron( (argc==1) ? FN : (strcmp("-",argv[1]) ? argv[1] : 0), &p);

  compVolumeIntegrals(&p);


  printf("\nT1 =   %+20.6f\n\n", T0);

  printf("Tx =   %+20.6f\n", T1[X]);
  printf("Ty =   %+20.6f\n", T1[Y]);
  printf("Tz =   %+20.6f\n\n", T1[Z]);

  printf("Txx =  %+20.6f\n", T2[X]);
  printf("Tyy =  %+20.6f\n", T2[Y]);
  printf("Tzz =  %+20.6f\n\n", T2[Z]);

  printf("Txy =  %+20.6f\n", TP[X]);
  printf("Tyz =  %+20.6f\n", TP[Y]);
  printf("Tzx =  %+20.6f\n\n", TP[Z]);

  density = 1.0;  /* assume unit density */

  mass = density * T0;

  /* compute center of mass */
  r[X] = T1[X] / T0;
  r[Y] = T1[Y] / T0;
  r[Z] = T1[Z] / T0;

  /* compute inertia tensor */
  J[X][X] = density * (T2[Y] + T2[Z]);
  J[Y][Y] = density * (T2[Z] + T2[X]);
  J[Z][Z] = density * (T2[X] + T2[Y]);
  J[X][Y] = J[Y][X] = - density * TP[X];
  J[Y][Z] = J[Z][Y] = - density * TP[Y];
  J[Z][X] = J[X][Z] = - density * TP[Z];

  /* translate inertia tensor to center of mass */
  J[X][X] -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
  J[Y][Y] -= mass * (r[Z]*r[Z] + r[X]*r[X]);
  J[Z][Z] -= mass * (r[X]*r[X] + r[Y]*r[Y]);
  J[X][Y] = J[Y][X] += mass * r[X] * r[Y];
  J[Y][Z] = J[Z][Y] += mass * r[Y] * r[Z];
  J[Z][X] = J[X][Z] += mass * r[Z] * r[X];

  printf("center of mass:  (%+12.6f,%+12.6f,%+12.6f)\n", r[X], r[Y], r[Z]);

  {
    double dpr = 180.0 / acos(-1.0);
    double x, y, z;
    double wlon, lat, radius;
    radius=sqrt(r[X]*r[X]+r[Y]*r[Y]+r[Z]*r[Z]);
    if ( radius == 0.0 ) { radius=1.0; }
    x=r[X]/radius;
    y=r[Y]/radius;
    z=r[Z]/radius;
    if ( x == 0.0 && y==0.0 ) { wlon=0.0; }
    else if ( x == 0.0 ) { wlon=y<0.0?90.0:-90.0; }
    else { wlon = atan2( -y, x) * dpr; }
    if ( wlon < 0.0) wlon += 360.0;
    lat = (z<-1.0) ? -999.0 : (z>1.0) ? 999.0 : (asin(z)*dpr);

    printf("  - lat,wlon,r:  (%+12.6f,%+12.6f,%+12.6f)\n\n", lat, wlon, radius);
  }

  printf("inertia tensor with origin at c.o.m. :\n");
  printf("%+15.6f  %+15.6g  %+15.6g\n", J[X][X], J[X][Y], J[X][Z]);
  printf("%+15.6g  %+15.6f  %+15.6g\n", J[Y][X], J[Y][Y], J[Y][Z]);
  printf("%+15.6g  %+15.6g  %+15.6f\n\n", J[Z][X], J[Z][Y], J[Z][Z]);

# define Ixx J[X][X]
# define Iyy J[Y][Y]
# define Izz J[Z][Z]
# define Ixy J[X][Y]
# define Iyz J[Y][Z]
# define Izx J[Z][X]

# define Iyx Ixy
# define Izy Iyz
# define Ixz Izx

  CC3 = -1.0;
  CC2 = Ixx + Iyy + Izz;
  CC1 = SQR(Ixy) + SQR(Iyz) + SQR(Izx)
        - (Ixx*Iyy + Iyy*Izz + Izz*Ixx);
  CC0 = (2.0*Ixy*Iyz*Izx) + (Ixx*Iyy*Izz)
        - (Izz*SQR(Ixy) + Ixx*SQR(Iyz) + Iyy*SQR(Izx));

  for (iEig=0; iEig<3; ++iEig) eigVals[iEig] = J[iEig][iEig];

  VolIntCubicRoots(CubicCoeffs, eigVals);

  printf( "Eigenvalues:\n    %35.25lg (%35.25lg)\n    %35.25lg (%35.25lg)\n    %35.25lg (%35.25lg)\n\n"
        , eigVals[0], EVALCUBIC(eigVals[0]), eigVals[1], EVALCUBIC(eigVals[1]), eigVals[2], EVALCUBIC(eigVals[2]));


  for (iEig=0; iEig<3; ++iEig)
  {
  double JI[3][3];
  int i, j, i1, j1, iDrop;
  int iMaxDet, i1MaxDet, jMaxDet, j1MaxDet;
  double det, maxDet;       /* maximum determinant magnitude */

    /*
    printf( "Matrix for eigenvalue %d = %lg:\n", iEig, eigVals[iEig]);
    /**/
    for ( j=0; j<3; ++j)
    {
      for ( i=0; i<3; ++i) JI[j][i] = J[j][i];
      JI[j][j] -= eigVals[iEig];
      /*
      for ( i=0; i<3; ++i) printf( " %19.12g", JI[j][i]);
      printf("\n");
      /**/
    }

    maxDet = 0.0;
    for ( i=0; i<2; ++i)
      for (j=0; j<2; ++j)
      {
        for ( i1=i+1; i1<3; ++i1) for (j1=j+1; j1<3; ++j1)
        {
          det = fabs(JI[j][i1]*JI[j1][i] - JI[j][i]*JI[j1][i1]);
          if ( det > maxDet)
          {
            maxDet = det;
            iMaxDet = i;
            i1MaxDet = i1;
            jMaxDet = j;
            j1MaxDet = j1;
          }
        }
      }

    if ( maxDet == 0.0)
    {
      printf( "Cannot solve for eigenvectors; probably X,Y,Z\n\n");
      return 0;
    }
    i = iMaxDet; i1 = i1MaxDet;
    j = jMaxDet; j1 = j1MaxDet;

    /* i=0; i1=1; j=0; j1=1;  /**/

    iDrop = 3 - (i + i1);

    det = JI[j][i]*JI[j1][i1] - JI[j][i1]*JI[j1][i];
    eigVecs[iEig][i1] = -(JI[j][i]*JI[j1][iDrop] - JI[j][iDrop]*JI[j1][i])/det;
    eigVecs[iEig][i] = (JI[j][i1]*JI[j1][iDrop] - JI[j][iDrop]*JI[j1][i1])/det;

    eigVecs[iEig][iDrop] = 1.0 / sqrt( 1.0 + SQR(eigVecs[iEig][i]) + SQR(eigVecs[iEig][i1]) );
    eigVecs[iEig][i1] *= eigVecs[iEig][iDrop];
    eigVecs[iEig][i] *= eigVecs[iEig][iDrop];
    /**/

    printf( "Eigenvec %d sol'n: ", iEig);
    for (j=0; j<3; ++j)
    {
    double result;
      result = 0.0;
      for (i=0; i<3; ++i) result += JI[j][i] * eigVecs[iEig][i];
      printf(" %+19.12g", result);
    }
    printf("\n");
  } /* for (iEig=0; iEig<3; ++iEig) */
  printf("\n");

  printf("Eigenvectors (in columns):\n");
  for ( iEig=0; iEig<3; ++iEig)
  {
    printf("%+15.12f %+15.12f %+15.12f\n"
          , eigVecs[X][iEig], eigVecs[Y][iEig], eigVecs[Z][iEig]);
  }
  printf("\n");

  printf("Eigenvectors in Lat & WLon (degrees):\n");
  for ( iEig=0; iEig<3; ++iEig)
  {
  double dpr = 180.0 / acos(-1.0);
  double y, x, z;
  double lat, wlon;
    x = eigVecs[iEig][X];
    y = eigVecs[iEig][Y];
    z = eigVecs[iEig][Z];

    lat = (z<-1.0) ? -999.0 : (z>1.0) ? 999.0 : asin(z);
    lat *= dpr;

    if ( x == 0.0 && y == 0.0)
    {
      wlon = 0.0;
    }
    else if ( x == 0.0)
    {
      wlon = (y<0 ? 90.0 : -90.0)/dpr;
    }
    else
    {
      wlon = atan2( -y, x);
    }
    wlon *= dpr;
    if ( wlon < 0.0) wlon += 360.0;

    printf("%+15.6f  %+15.6fW   (Lat, WLon in degrees)\n", lat, wlon);
  }
  printf("\n");

  return 0;
}
