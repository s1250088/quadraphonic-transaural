#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>

#define PDMAC_VERSION

#include "m_pd.h"
#include "fftw3.h"
#include "sqlite3.h"

#define DEFAULT_N_POINTS 256
#define MAX_N_POINTS 3000
#define MAX_BLOCKSIZE 8192
#define MAX_SAMPLES 4096
#define F_OK    0
#define DBFILENAME "/HRIR@44100.db"
#define MyPI 3.14159265358979323846


static t_class *quadraTrans_tilde_class;

typedef struct _quadraTrans_tilde{
    t_object x_obj;
    t_sample f_quadraTrans;
    t_sample f;
    
    t_float aziL;
    
    t_float sr;
    
    t_inlet *x_in2;
    t_inlet *x_in3;
    
    t_outlet *x_outL;
    t_outlet *x_outR;
    t_outlet *x_outSL;
    t_outlet *x_outSR;
    
    t_float previousImpulse[2][MAX_N_POINTS];
    t_float currentImpulse[9][2][MAX_N_POINTS];
    int currentRow;
    t_float convBufferC[MAX_N_POINTS];
    t_float convBufferI[MAX_N_POINTS];
    
    //fftw
    int fftsize;
    fftwf_plan fftplanL;
    float *fftinL;
    fftwf_complex *fftoutL;
    fftwf_plan fftplanR;
    float *fftinR;
    fftwf_complex *fftoutR;
    
    //invfftw
    fftwf_plan fftplanInvL;
    fftwf_complex *fftinInvL;
    float *fftoutInvL;
    fftwf_plan fftplanInvR;
    fftwf_complex *fftinInvR;
    float *fftoutInvR;
    
    //hrirfft
    fftwf_plan fftplanHrirL;
    float *fftinHrirL;
    fftwf_complex *fftoutHrirL;
    fftwf_plan fftplanHrirR;
    float *fftinHrirR;
    fftwf_complex *fftoutHrirR;
    
    sqlite3 *db;
    t_int connected;
    t_int rc;
    
    sqlite3 *hrtfDb;
    t_int hrtfConnected;
    t_int hrtfRc;
    
}t_quadraTrans_tilde;

t_int *quadraTrans_tilde_perform(t_int *w){
    t_quadraTrans_tilde *x = (t_quadraTrans_tilde *)(w[1]);
    t_sample  *inXL   =  (t_sample *)(w[2]);
    t_sample  *inXR   =  (t_sample *)(w[3]);
    t_sample  *outSR =  (t_sample *)(w[4]);
    t_sample  *outSL =  (t_sample *)(w[5]);
    t_sample  *outR  =  (t_sample *)(w[6]);
    t_sample  *outL  =  (t_sample *)(w[7]);
    int           n  =         (int)(w[8]);
    
    int az = (int)(x->aziL);
    
    float tmp;
    for (int i = 0; i < n + 512 -1; i++) {
        if (i < n) {
            tmp = inXL[i];
            x->fftinL[i] = ((inXR[i] + inXL[i]))/2;
            x->fftinR[i] = ((tmp - inXR[i]) )/2;
        } else {
            x->fftinL[i] = 0;
            x->fftinR[i] = 0;
        }
        
        if (i < 512) {
            //x->fftinHrirL[i] = x->currentImpulse[x->currentRow][0][i];
            //x->fftinHrirR[i] = x->currentImpulse[x->currentRow][1][i];
            x->fftinHrirL[i] = 1 / (x->currentImpulse[x->currentRow][0][i] + x->currentImpulse[x->currentRow][1][i]);
            x->fftinHrirR[i] = 1 / (x->currentImpulse[x->currentRow][0][i] - x->currentImpulse[x->currentRow][1][i]);
            if (x->fftinHrirL[i] <= 0.000001 && x->fftinHrirR[i] <= 0.000001) {
                //x->fftinHrirL[i] = 0.00001;
                //x->fftinHrirR[i] = 0.00001;
            }
        } else {
            x->fftinHrirL[i] = 0;
            x->fftinHrirL[i] = 0;
        }
    }
    
    //FFT
    fftwf_execute(x->fftplanL);
    fftwf_execute(x->fftplanR);
    fftwf_execute(x->fftplanHrirL);
    fftwf_execute(x->fftplanHrirR);
    
    //multiple
    //fftw_complex a[DEFAULT_N_POINTS];
    //fftw_complex b[DEFAULT_N_POINTS];
    
    for (int i = 0; i < (n + 512 - 1)/2+1; i++) {
        
        x->fftinInvL[i][0] = x->fftoutL[i][0] * x->fftoutHrirL[i][0] - x->fftoutL[i][1] * x->fftoutHrirL[i][1];
        x->fftinInvL[i][1] = x->fftoutL[i][0] * x->fftoutHrirL[i][1] + x->fftoutL[i][1] * x->fftoutHrirL[i][0];
        
        x->fftinInvR[i][0] = x->fftoutR[i][0] * x->fftoutHrirR[i][0] - x->fftoutR[i][1] * x->fftoutHrirR[i][1];
        x->fftinInvR[i][1] = x->fftoutR[i][0] * x->fftoutHrirR[i][1] + x->fftoutR[i][1] * x->fftoutHrirR[i][0];
    }
    
    // inv FFT
    fftwf_execute(x->fftplanInvL);
    fftwf_execute(x->fftplanInvR);
    
    // Outputs
    for (int i = 0; i < n; i++) {
        tmp = x->fftoutInvL[i];
        x->fftoutInvL[i] += x->fftoutInvR[i];
        x->fftoutInvR[i] = tmp - x->fftoutInvR[i];
        
        if(az<360) az += 360;
        else if(az>360) az = az%360;
        
        if(315<az || az<45){
            *outL++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outR++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(az==45){
            *outL++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outSR++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(45<az || az<135){
            *outR++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outSR++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(az==135){
            *outR++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outSL++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(135<az || az<225){
            *outSR++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outSL++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(az==225){
            *outSR++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outL++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(225<az || az<315){
            *outSL++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outL++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
        else if(az==315){
            *outSL++ = x->fftoutInvL[i] / (x->fftsize*100000);
            *outR++ = x->fftoutInvR[i] / (x->fftsize * 100000);
        }
    }
    
    return (w+9);
}

void quadraTrans_tilde_dsp(t_quadraTrans_tilde *x, t_signal **sp){
    dsp_add(quadraTrans_tilde_perform, 8, x,
            sp[0]->s_vec, sp[1]->s_vec,
            sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
            sp[0]->s_n);
    
    x->fftsize = (sp[0]->s_n+ 512 -1); //to have a blocksize fft output
    
    //set L
    x->fftinL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanL = fftwf_plan_dft_r2c_1d(x->fftsize, x->fftinL, x->fftoutL, FFTW_MEASURE);
    //set R
    x->fftinR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinR, x->fftoutR, FFTW_MEASURE);
    //set invL
    x->fftinInvL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftoutInvL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftplanInvL = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvL, x->fftoutInvL, FFTW_MEASURE);
    //set invR
    x->fftinInvR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftoutInvR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftplanInvR = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvR, x->fftoutInvR, FFTW_MEASURE);
    //set Hrir
    x->fftinHrirL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutHrirL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanHrirL = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirL, x->fftoutHrirL, FFTW_MEASURE);
    x->fftinHrirR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutHrirR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanHrirR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirR, x->fftoutHrirR, FFTW_MEASURE);
}

void quadraTrans_tilde_free(t_quadraTrans_tilde *x){
    //inlet
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    
    //outlet
    outlet_free(x->x_outL);
    outlet_free(x->x_outR);
    outlet_free(x->x_outSL);
    outlet_free(x->x_outSR);
    
    //fftw
    fftwf_free(x->fftoutL);
    fftwf_free(x->fftinL);
    fftwf_destroy_plan(x->fftplanL);
    fftwf_free(x->fftoutR);
    fftwf_free(x->fftinR);
    fftwf_destroy_plan(x->fftplanR);
    fftwf_free(x->fftoutInvL);
    fftwf_free(x->fftinInvL);
    fftwf_destroy_plan(x->fftplanInvL);
    fftwf_free(x->fftoutInvR);
    fftwf_free(x->fftinInvR);
    fftwf_destroy_plan(x->fftplanInvR);
    fftwf_free(x->fftoutHrirL);
    fftwf_free(x->fftinHrirL);
    fftwf_destroy_plan(x->fftplanHrirL);
    fftwf_free(x->fftoutHrirR);
    fftwf_free(x->fftinHrirR);
    fftwf_destroy_plan(x->fftplanHrirR);
   
}

void *quadraTrans_tilde_new(t_floatarg f){
    t_quadraTrans_tilde *x = (t_quadraTrans_tilde *)pd_new(quadraTrans_tilde_class);
    
    x->f_quadraTrans = f;
    
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->aziL);
    
    x->x_outL  = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outR  = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outSL = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outSR = outlet_new(&x->x_obj, gensym("signal"));
    
    return (void *)x;
}

void quadraTrans_tilde_setup(void) {
    quadraTrans_tilde_class = class_new(gensym("quadraTrans~"),
                                        (t_newmethod)quadraTrans_tilde_new,
                                        (t_method)quadraTrans_tilde_free,
                                        sizeof(t_quadraTrans_tilde),
                                        CLASS_DEFAULT,
                                        A_DEFFLOAT, 0);
    
    class_addmethod(quadraTrans_tilde_class, (t_method)quadraTrans_tilde_dsp, gensym("dsp"), 0);
    CLASS_MAINSIGNALIN(quadraTrans_tilde_class, t_quadraTrans_tilde, f);
}