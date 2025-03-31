#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>

#include "cpgplot.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include "gsl/gsl_sort.h"
#include "gsl/gsl_statistics.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multifit_nlinear.h"

#define PI 3.141592653589793238462643383279
#define CLINE 299792458

void gaussian_fit(float *data_t,float *data,int n,float *data_fit,double *a,double *b,double *c);
double gaussian(const double a, const double b, const double c, const double t);
int func_f (const gsl_vector * x, void *params, gsl_vector * f);
int func_df (const gsl_vector * x, void *params, gsl_matrix * J);
int func_fvv (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv);
void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
void solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,gsl_multifit_nlinear_parameters *params);

  struct data
  {
    double *t;
    double *y;
    size_t n;
  };

void help();
int main(int argc, char* argv[])
{
    int N=0;
    int i;
    char filename[1000];
    FILE *fp;

    if(argc < 2) help();
    for(i=0;i<argc;i++)
    {
        if (strcmp(argv[i],"-h")==0) help();
    }
    strcpy(filename,argv[argc-1]);//获取文件名
    printf("Reading from file: %s\n",filename);

    
    fp=fopen(filename,"r");
    while(!feof(fp))                  //获取list中的p0数
    {
        if(fgetc(fp)=='\n')
        {
            N++;
        }
    }
    printf("P0 num: %d\n",N);
    fclose(fp);

    double  p0array0[N];
    double  accarray0[N];
    
    float  p0array[N];
    float  accarray[N];
    float  accarray2[N];
    float  accarrayf[N];

    fp=fopen(filename,"r");
    for ( i = 0; i < N; i++)
    {
        fscanf(fp,"%lf %lf",&p0array0[i],&accarray0[i]);
        p0array[i]=p0array0[i]/1000;                        //ms convert to s
        accarray[i]=accarray0[i];
        accarray2[i]=accarray0[i]*accarray0[i];
    }
    fclose(fp);


    // float minp0=p0array[0],maxp0=p0array[0];
    // float maxacc=fabs(accarray0[0]);
    // float minacc;

    // for ( i = 1; i < N; i++)
    // {
    //     if(minp0>p0array[i]) minp0=p0array[i];
    //     if(maxp0<p0array[i]) maxp0=p0array[i];
    //     if(maxacc<fabs(accarray0[i])) maxacc=fabs(accarray0[i]);
    // }
    // minacc=0-maxacc;
    
    double A,B,C;
    gaussian_fit(p0array,accarray2,N,accarrayf,&A,&B,&C);
    float P0=-1.0*B/(2*A);  
    float P1=sqrt(fabs(C/A-B*B/4.0/A/A));
    float A1=P1*sqrt(-1.0*A);
    float PB=2*PI*CLINE/P0/sqrt(-1.0*A);
    float x=P1*P1/P0/P0*CLINE/A1;
    // PB=P1/P0*2*PI*CLINE/A1; // same as PB=2*PI*CLINE/P0/sqrt(-1.0*A);
    
    printf("P0 = %f us\n",P0*1000000);
    printf("P1 = %f us\n",P1*1000000);
    printf("A1 = %f m*s^-2\n",A1);
    printf("PB = %.15f min\n",PB/60);
    printf("x=ap*sin(i)/c = %f s\n",x);

    for(i=0;i<N;i++)
    {
        p0array[i]*=1000000;        // s convert to us
    }

    int pointP=501;
    float *largA=(float *)malloc(sizeof(float)*pointP);
    float *largP=(float *)malloc(sizeof(float)*pointP);
    largP[0]=(P0-P1)*1000000;
    largA[0]=0;
    largP[(pointP-1)/2]=(P0+P1)*1000000;
    largA[(pointP-1)/2]=0;
    for(i=1;i<(pointP-1)/2;i++)
    {
        // float thisP0=P0-P1+i*4.0*P1/(pointP-1);
        largA[i] = sqrt(gaussian(A,B,C,P0-P1+i*4.0*P1/(pointP-1)));
        // largA[i] = sqrt(A*thisP0*thisP0 + B*thisP0 + C);
        largA[pointP-1-i]=0.0-largA[i];

        largP[i] = (P0-P1+i*4*P1/(pointP-1))*1000000;
        largP[pointP-1-i] = largP[i];
    }
    largP[pointP-1]=(P0-P1)*1000000;
    largA[pointP-1]=0;

    // for(i=0;i<pointP;i++)
    // {
    //     printf("A=%f P=%f\n",largA[i],largP[i]);
    // }

    float OP[1];
    float OA[1];
    char txt0[1000];
    char txtx[1000];
    char txty[1000];
    OP[0]=P0*1000000;
    OA[0]=0;
    
    sprintf(txtx,"%f",P1*1000000);
    sprintf(txty,"%.2f",A1);

    char outtxt[1000];
    sprintf(outtxt,"%s.ps/cps",filename);
    // cpgbeg(0,"/xs",1,1);
    cpgbeg(0,outtxt,1,1);
    cpgsvp(0.1,0.9,0.1,0.9);
    cpgslw(2);
    // cpgswin(P0*1000000-P1*1000000*1.2,P0*1000000+P1*1000000*1.2,minacc-(maxacc-minacc)*0.1,maxacc+(maxacc-minacc)*0.1);
    cpgswin(P0*1000000-P1*1000000*1.2,P0*1000000+P1*1000000*1.2,0-A1*1.2,0+A1*1.2);
    cpgpt(N,p0array,accarray,2);
    cpgpt(1,OP,OA,9);
    sprintf(txt0,"P0: %f us",P0*1000000);
    cpgptxt(P0*1000000,0-A1*0.12,0,0.5,txt0);
    sprintf(txt0,"Pb: %.2f min",PB/60);
    cpgptxt(P0*1000000,0-A1*0.22,0,0.5,txt0);
    sprintf(txt0,"x: %f s",x);
    cpgptxt(P0*1000000,0-A1*0.32,0,0.5,txt0);


    cpgarro(OP[0],0.0,OP[0]+P1*1000000,0.0);
    cpgptxt(P0*1000000+P1*1000000*0.5,0+A1*0.05,0,0.5,txtx);
    cpgarro(OP[0],0.0,OP[0],A1);
    cpgptxt(P0*1000000-P1*1000000*0.05,0+A1*0.5,90,0.5,txty);
    
    cpgsci(4);
    //cpgline(pointP-4,largP+2,largA+2);
    cpgline(pointP,largP,largA);
    cpgsci(1);
    cpgbox("bcnst",0,0,"bcnst",0,0);
    cpglab("P0 (us)","Acceleration (m s\\u-2\\d)","");
    cpgend();



    free(largA);
    free(largP);
    return 0;
}

void gaussian_fit(float *data_t,float *data,int n,float *data_fit,double *a,double *b,double *c)
{
    const size_t p = 3;
    const gsl_rng_type * T = gsl_rng_default;
    gsl_vector *f = gsl_vector_alloc(n);
    gsl_vector *x = gsl_vector_alloc(p);
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
    struct data fit_data;
    gsl_rng * r;
    size_t i;
    gsl_rng_env_setup ();
    r = gsl_rng_alloc (T);
    fit_data.t = malloc(n * sizeof(double));
    fit_data.y = malloc(n * sizeof(double));
    fit_data.n = n;
    /* generate synthetic data with noise */
    for (i = 0; i < n; ++i)
    {
        //double t = (double)i / (double) n;
        //double y0 = gaussian(a, b, c, t);
        //double dy = gsl_ran_gaussian (r, 0.1 * y0);
        fit_data.t[i] = data_t[i];
        fit_data.y[i] = data[i];
    }
    /* define function to be minimized */
    fdf.f = func_f;
    fdf.df = func_df;
    fdf.fvv = func_fvv;
    // fdf.fvv = NULL;
    fdf.n = n;
    fdf.p = p;
    fdf.params = &fit_data;
    /* starting point */
    gsl_vector_set(x, 0, 1.0);
    gsl_vector_set(x, 1, 0.0);
    gsl_vector_set(x, 2, 1.0);
    fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    solve_system(x, &fdf, &fdf_params);
    /* print data and model */
    {
    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);
    for (i = 0; i < n; ++i)
    {
        double ti = fit_data.t[i];
        //double yi = fit_data.y[i];
        double fi = gaussian(A, B, C, ti);
        //printf("%f %f %f\n", ti, yi, fi);
        data_fit[i]=fi;
    }
    *a=A;
    *b=B;
    *c=C;
    }
    gsl_vector_free(f);
    gsl_vector_free(x);
    gsl_rng_free(r);
}

// /* model function: a * exp( -1/2 * [ (t - b) / c ]^2 ) */
/*y=a*x^2 + b*x + c*/
double gaussian(const double a, const double b, const double c, const double t)
{
    // const double z = (t - b) / c;
    // return (a * exp(-0.5 * z * z));

    return (a*t*t+b*t+c);
}
int func_f (const gsl_vector * x, void *params, gsl_vector * f)
{
    struct data *d = (struct data *) params;
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double c = gsl_vector_get(x, 2);
    size_t i;
    for (i = 0; i < d->n; ++i)
    {
        double ti = d->t[i];
        double yi = d->y[i];
        double y = gaussian(a, b, c, ti);
        gsl_vector_set(f, i, yi - y);
    }
    return GSL_SUCCESS;
}
int func_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
    struct data *d = (struct data *) params;
    // double a = gsl_vector_get(x, 0);
    // double b = gsl_vector_get(x, 1);
    // double c = gsl_vector_get(x, 2);
    size_t i;
    for (i = 0; i < d->n; ++i)
    {
        // double ti = d->t[i];
        // double zi = (ti - b) / c;
        // double ei = exp(-0.5 * zi * zi);
        gsl_matrix_set(J, i, 0, -1.0*d->t[i]*d->t[i]);
        gsl_matrix_set(J, i, 1, -1.0*d->t[i]);
        gsl_matrix_set(J, i, 2, -1.0);
    }
    return GSL_SUCCESS;
}
int func_fvv (const gsl_vector * x, const gsl_vector * v, void *params, gsl_vector * fvv)
{
    struct data *d = (struct data *) params;
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double c = gsl_vector_get(x, 2);
    double va = gsl_vector_get(v, 0);
    double vb = gsl_vector_get(v, 1);
    double vc = gsl_vector_get(v, 2);
    size_t i;
    for (i = 0; i < d->n; ++i)
    {
        double ti = d->t[i];
        double zi = (ti - b) / c;
        double ei = exp(-0.5 * zi * zi);
        double Dab = -zi * ei / c;
        double Dac = -zi * zi * ei / c;
        double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
        double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
        double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
        double sum;
        sum = 2.0 * va * vb * Dab +
        2.0 * va * vc * Dac +
        vb * vb * Dbb +
        2.0 * vb * vc * Dbc +
        vc * vc * Dcc;
        sum=0;
        gsl_vector_set(fvv, i, sum);
    }
    return GSL_SUCCESS;
}
void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    //gsl_vector *f = gsl_multifit_nlinear_residual(w);
    //gsl_vector *x = gsl_multifit_nlinear_position(w);
    //double avratio = gsl_multifit_nlinear_avratio(w);
    //double rcond;
    (void) params; /* not used */
    /* compute reciprocal condition number of J(x) */
    //gsl_multifit_nlinear_rcond(&rcond, w);
    /*fprintf(stderr, "iter %2zu: a = %.4f, b = %.4f, c = %.4f, |a|/|v| = %.4f cond(J) =
    ,!%8.4f, |f(x)| = %.4f\n",iter,gsl_vector_get(x, 0),gsl_vector_get(x, 1),gsl_vector_get(x, 2),avratio,1.0 / rcond,gsl_blas_dnrm2(f));
    */
}
void solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,gsl_multifit_nlinear_parameters *params)
{
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    const size_t max_iter = 1000;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    const double ftol = 1.0e-8;
    const size_t n = fdf->n;
    const size_t p = fdf->p;
    gsl_multifit_nlinear_workspace *work =
    gsl_multifit_nlinear_alloc(T, params, n, p);
    gsl_vector * f = gsl_multifit_nlinear_residual(work);
    gsl_vector * y = gsl_multifit_nlinear_position(work);
    int info;
    double chisq0, chisq, rcond;
    /* initialize solver */
    gsl_multifit_nlinear_init(x, fdf, work);
    /* store initial cost */
    gsl_blas_ddot(f, f, &chisq0);
    /* iterate until convergence */
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
    callback, NULL, &info, work);
    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);
    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);
    gsl_vector_memcpy(x, y);
    /* print summary */
    //fprintf(stderr, "NITER = %zu\n", gsl_multifit_nlinear_niter(work));
    //fprintf(stderr, "NFEV = %zu\n", fdf->nevalf);
    //fprintf(stderr, "NJEV = %zu\n", fdf->nevaldf);
    //fprintf(stderr, "NAEV = %zu\n", fdf->nevalfvv);
    //fprintf(stderr, "initial cost = %.12e\n", chisq0);
    //fprintf(stderr, "final cost = %.12e\n", chisq);
    //fprintf(stderr, "final x = (%.12e, %.12e, %12e)\n",
    //gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2));
    //fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);
    gsl_multifit_nlinear_free(work);
}

void help()
{
    printf(" This program is build for fit P0-Acc of binary pulsar data\n");
    printf(" Usage: \n");
    printf("        p0afit   filename\n");
    printf(" The input file should be a two col file that include P0(ms) for the 1st col\n and Acc (m s^-2) for the 2nd col\n");
    exit(0);
}
