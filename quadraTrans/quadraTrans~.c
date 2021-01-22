#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>

//#define PDMAC_VERSION

#include "m_pd.h"
#include "fftw3.h"
#include "sqlite3.h"

#define DEFAULT_N_POINTS 256
#define MAX_N_POINTS 3000
#define MAX_BLOCKSIZE 8192
#define MAX_SAMPLES 4096
#define F_OK    0
#define DBFILENAME "./HRIR@44100.db"
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
    
    t_int nPts;
    t_int nPtsInDb;
    t_float azimuth;
    t_float elevation;
    t_float distance;
    
    t_float z;
    t_float g;
    t_int   m;
    
    t_float azi;
    t_float ele;
    t_float dis;
    t_float rs[2];
    t_float es[2];
    t_float as[2];
    t_float on[3];
    t_float set[9][3];
    t_float dist[9];
    t_int n_meas;
    
    t_float prevAzi;
    t_float prevEle;
    t_float prevRng;
    double *window;
    float hrtf[MAX_N_POINTS];
    
    t_float OUT[2][MAX_N_POINTS];
    
    char *zErrMsg;
    char *zHrtfErrMsg;
    
    t_float previousImpulse[2][MAX_N_POINTS];
    t_float currentImpulse[9][2][MAX_N_POINTS];
    int currentRow;
    t_float convBufferC[MAX_N_POINTS];
    t_float convBufferI[MAX_N_POINTS];
    
    t_float low_filter_hz;
    
    char path[2000];
    t_float crossCoef[MAX_BLOCKSIZE];
    t_int bufferPin;
    
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

#include "quadraTransFunctions.h"

static void makewindow(double *w, int n) {
    int i;
    double xshift = n / 2.0;
    double x;
    for (i = 0; i < n; i++) {
        x = (i - xshift) / xshift;
        w[i] = 0.5 * (1 + cos(MyPI * x));
    }
}

t_int *quadraTrans_tilde_perform(t_int *w){
    t_quadraTrans_tilde *x = (t_quadraTrans_tilde *)(w[1]);
    t_sample  *inXL   =  (t_sample *)(w[2]);
    t_sample  *inXR   =  (t_sample *)(w[3]);
    t_sample  *outSR =  (t_sample *)(w[4]);
    t_sample  *outSL =  (t_sample *)(w[5]);
    t_sample  *outR  =  (t_sample *)(w[6]);
    t_sample  *outL  =  (t_sample *)(w[7]);
    int   blocksize  =         (int)(w[8]);
    
    int az = (int)(x->aziL);
    
    if (x->connected == 1) {
        x->prevAzi = x->azi;
        x->prevEle = x->ele;
        x->prevRng = x->dis;
        
        x->azi = roundf(x->azimuth);
        x->ele = roundf(x->elevation);
        x->dis = roundf(x->distance);
        if ((x->prevAzi == x->azi) && (x->prevEle == x->ele) && (x->prevRng == x->dis)) {
         } else {
         rangeDis(x);
         rangeEle(x);
         rangeAzi(x);
         formSet(x);
         testCoplan_lin(x);
         }
        
#pragma region main
        for (int i = 0; i < blocksize; i++) {
            if (i < blocksize) {
                x->fftinL[i] = (inXL[i] + inXR[i])/2;
                x->fftinR[i] = (inXL[i] - inXR[i])/2;
            } /*else {
                x->fftinL[i] = 0;
                x->fftinR[i] = 0;
            }*/
            
            fftwf_execute(x->fftplanL);
            fftwf_execute(x->fftplanR);
            
            if (i < blocksize) {
                x->fftinHrirL[i] = 1 / (x->currentImpulse[x->currentRow][0][i] + x->currentImpulse[x->currentRow][1][i]);
                x->fftinHrirR[i] = 1 / (x->currentImpulse[x->currentRow][0][i] - x->currentImpulse[x->currentRow][1][i]);
            } /*else {
                x->fftinHrirL[i] = 0;
                x->fftinHrirL[i] = 0;
            }*/
            fftwf_execute(x->fftplanHrirL);
            fftwf_execute(x->fftplanHrirR);
        }
        
        for (int i = 0; i < blocksize; i++) {
            
            x->fftinInvL[i][0] = x->fftoutL[i][0] * x->fftoutHrirL[i][0] - x->fftoutL[i][1] * x->fftoutHrirL[i][1];
            x->fftinInvL[i][1] = x->fftoutL[i][0] * x->fftoutHrirL[i][1] + x->fftoutL[i][1] * x->fftoutHrirL[i][0];
            
            x->fftinInvR[i][0] = x->fftoutR[i][0] * x->fftoutHrirR[i][0] - x->fftoutR[i][1] * x->fftoutHrirR[i][1];
            x->fftinInvR[i][1] = x->fftoutR[i][0] * x->fftoutHrirR[i][1] + x->fftoutR[i][1] * x->fftoutHrirR[i][0];
        }
        
        fftwf_execute(x->fftplanInvL);
        fftwf_execute(x->fftplanInvR);
        
        // Outputs
        for (int i = 0; i < blocksize; i++) {
            float tmp;
            tmp = x->fftoutInvL[i];
            x->fftoutInvL[i] += x->fftoutInvR[i];
            x->fftoutInvR[i] = tmp - x->fftoutInvR[i];
            
            if(az<360) az += 360;
            else if(az>360) az = az%360;
            
            if(315<az || az<45){
                *outL++ = x->fftoutInvL[i] / x->fftsize;
                *outR++ = x->fftoutInvR[i] / x->fftsize;
            }
            else if(az==45){
                *outL++  = x->fftoutInvL[i] / x->fftsize;
                *outSR++ = x->fftoutInvR[i] / x->fftsize;
            }
            else if(45<az || az<135){
                *outR++  = x->fftoutInvL[i] / x->fftsize;
                *outSR++ = x->fftoutInvR[i] / x->fftsize;
            }
            else if(az==135){
                *outR++  = x->fftoutInvL[i] / x->fftsize;
                *outSL++ = x->fftoutInvR[i] / x->fftsize;
            }
            else if(135<az || az<225){
                *outSR++ = x->fftoutInvL[i] / x->fftsize;
                *outSL++ = x->fftoutInvR[i] / x->fftsize;
            }
            else if(az==225){
                *outSR++ = x->fftoutInvL[i] / x->fftsize;
                *outL++  = x->fftoutInvR[i] / x->fftsize;
            }
            else if(225<az || az<315){
                *outSL++ = x->fftoutInvL[i] / x->fftsize;
                *outL++  = x->fftoutInvR[i] / x->fftsize;
            }
            else if(az==315){
                *outSL++ = x->fftoutInvL[i] / x->fftsize;
                *outR++  = x->fftoutInvR[i] / x->fftsize;
            }
        }
#pragma endregion
    }
    return (w+9);
}

void quadraTrans_tilde_dsp(t_quadraTrans_tilde *x, t_signal **sp){
    dsp_add(quadraTrans_tilde_perform, 8, x,
            sp[0]->s_vec, sp[1]->s_vec,
            sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
            sp[0]->s_n);
    
#pragma region fftw_set
    x->fftsize = (sp[0]->s_n);
    
    //set L
    x->fftinL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanL = fftwf_plan_dft_r2c_1d(x->fftsize, x->fftinL, x->fftoutL, FFTW_MEASURE);
    //set R
    x->fftinR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinR, x->fftoutR, FFTW_MEASURE);
    //set HrirL
    x->fftinHrirL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutHrirL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanHrirL = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirL, x->fftoutHrirL, FFTW_MEASURE);
    //set HrirR
    x->fftinHrirR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftoutHrirR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftplanHrirR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirR, x->fftoutHrirR, FFTW_MEASURE);
    //set invL
    x->fftinInvL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftoutInvL = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftplanInvL = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvL, x->fftoutInvL, FFTW_MEASURE);
    //set invR
    x->fftinInvR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->fftoutInvR = (float *) fftwf_malloc(sizeof(float) * x->fftsize);
    x->fftplanInvR = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvR, x->fftoutInvR, FFTW_MEASURE);
    
#pragma endregion fftw_set
    
#pragma region db
    x->window = getbytes(sizeof(double)* x->fftsize);
    makewindow(x->window, sp[0]->s_n);
    int power, base;
    x->sr = sp[0]->s_sr;
    char file[2000] = "";
    char str[8] = "";
    
    strcpy(file,  x->path);
    strcat(file,  "./HRIR@");
    sprintf(str,  "%d", (int)x->sr);
    strcat(file,  str);
    strcat(file,  ".db");
    if (access(file, F_OK) == -1) {
        post("Warning: The database %s doesn't exist,\ntrying to use %s database", file, DBFILENAME);
        
        strcpy(file, x->path);
        strcat(file,  "./HRIR@");
        sprintf(str,  "%d", 44100);
        strcat(file,  str);
        strcat(file,  ".db");
        
        if (access(file, F_OK) == -1)
            post("Error: The database %s doesn't exist either, trying downloading from \n http://arts.u-aizu.ac.jp/research/hrir-with-distance-control/");
        else post("Warning: Using %s at a sampling freq. of %d, results are not correct.", file, (int)x->sr);
    } else {
        if (x->connected == 1) {
            sqlite3_close(x->db);
            x->connected = 0;
        }
        x->rc = sqlite3_open(file, &(x->db));
        if (x->rc) {
            post("Can't open database: %s : %s\n", file, sqlite3_errmsg(x->db));
            sqlite3_close(x->db);
        } else {
            post("%s database opened", file);
            x->connected = 1;
        }
        
        switch ((int)(x->sr)) {
            case 8000: // number of taps at 8000 = 125
                x->nPtsInDb = 250;
                break;
            case 16000: // number of taps at 16000 = 250
                x->nPtsInDb = 500;
                break;
            case 22050: //number of taps at 22050 = 344
                x->nPtsInDb = 688;
                break;
            case 44100: // number of taps at 44100 = 689
                x->nPtsInDb = 1378;
                break;
            case 48000: // number of taps at 48000 = 750
                x->nPtsInDb = 1500;
                break;
            case 65536: // number of taps at 65536 = 1024
                x->nPtsInDb = 2048;
                break;
            case 96000: // number of taps at 96000 = 1500
                x->nPtsInDb = 3000;
                break;
            case 192000: // number of taps at 192000 = 3000
                x->nPtsInDb = 6000;
                break;
            default:
                post("WARNING: Define the right amount of hrir taps for the  sampling rate: %0.0f Hz.\nConsider downloading a suitable database for this sampling rate.\n", x->sr);
                break;
        }
        
        base = x->nPtsInDb / 2;
        power = 0;
        
        while (base > 1) {
            base = base >> 1;
            power += 1;
        }
        base = 1;
        while (power-- > 0) base = base << 1;
        
        if (x->nPts > x->nPtsInDb / 2) x->nPts = base;
        
        if (!isPowerOfTwo(x->nPts)) {
            x->nPts = base;
            post("Number of tabs must be a power of 2\n", base, x->nPts);
        }
        post("hrir~: Max. blocksize: 8192, Max. taps %d, currently using %d \n", base, x->nPts);
    }
#pragma endregion db
}

void quadraTrans_tilde_free(t_quadraTrans_tilde *x){
    freebytes(x->window, sizeof(double)* x->fftsize);
    
    //inlet
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    
    //outlet
    outlet_free(x->x_outL);
    outlet_free(x->x_outR);
    outlet_free(x->x_outSL);
    outlet_free(x->x_outSR);
    
    if (x->connected) {
        sqlite3_close(x->db);
        post("measurement database closed");
        x->connected = 0;
    }
    
#pragma region fftw_free
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
#pragma endregion fftw_free
    
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
    
#pragma region init
    x->azi = 45;
    x->ele = 0;
    x->dis = 140;
    x->nPtsInDb = 300;
    x->elevation = 0;
    x->azimuth = 0;
    x->distance = 20;
    
    x->nPts = 512;
    
    for (int i = 0; i < MAX_N_POINTS; i++) {
        x->convBufferC[i] = 0;
        x->convBufferI[i] = 0;
        x->previousImpulse[0][i] = 0;
        x->previousImpulse[1][i] = 0;
        
        x->OUT[0][i] = 0;
        x->OUT[1][i] = 0;
        
        x->crossCoef[i] = 1.0 * i / (MAX_BLOCKSIZE - 1.0);
        x->crossCoef[i] = cos((MyPI / 2) * ((float)i / (MAX_BLOCKSIZE - 1.0)));
    }
    
    strcat(x->path, canvas_getdir(canvas_getcurrent())->s_name);
    x->connected = 0;
#pragma endregion init
    
    return (void *)x;
}

void quadraTrans_tilde_setup(void) {
    quadraTrans_tilde_class = class_new(gensym("quadraTrans~"),
                                        (t_newmethod)quadraTrans_tilde_new,
                                        (t_method)quadraTrans_tilde_free,
                                        sizeof(t_quadraTrans_tilde),
                                        CLASS_DEFAULT,
                                        A_DEFSYMBOL, A_DEFSYMBOL, 0);
    
    class_addmethod(quadraTrans_tilde_class, (t_method)quadraTrans_tilde_dsp, gensym("dsp"), 0);
    CLASS_MAINSIGNALIN(quadraTrans_tilde_class, t_quadraTrans_tilde, f);
}
