#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>

#define PDMAC_VERSION
#ifdef WIN_VERSION
#include <io.h>
#endif

#include "m_pd.h"
#include "fftw3.h"
#include "sqlite3.h"

#define DEFAULT_N_POINTS 256
#define MAX_N_POINTS 3000
#define MAX_BLOCKSIZE 8192
#define DBFILENAME "/HRIR@44100.db"

#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)
#define radToDeg(angleInRadians) ((angleInRadians) * 180.0 / M_PI)

static t_class *quadraTrans_tilde_class;

typedef struct _quadraTrans_tilde{
    t_object x_obj;
    t_sample f_quadraTrans;
    t_float f;
    
    t_float aziL; //listener azimuth
    
    int sr; //sampling rate
    
    t_inlet *x_in2; //inlet
    t_inlet *x_in3; //inlet
    
    t_outlet *x_outL; //outlet
    t_outlet *x_outR; //outlet
    t_outlet *x_outSL; //outlet
    t_outlet *x_outSR; //outlet
    
    int convsize;
    int nbins;
    
    int nPts; //No. of taps used for the convolution
    int nPtsInDb; //No. of taps existing in the database (sampling rate dependent)
    
    t_float currentImpulse[2][MAX_N_POINTS];
    t_float speackersFilters[4][2][MAX_N_POINTS]; //0->L; 1->R; 2->SL; 3->SR    
    t_sample buffer[MAX_BLOCKSIZE*MAX_N_POINTS][2];
    char path[2000];
    
    sqlite3 *db;
    t_int connected;
    t_int rc;
    char *zErrMsg;
    
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
    
}t_quadraTrans_tilde;

#include "quadraTransFunctions.h"

t_int *quadraTrans_tilde_perform(t_int *w){
    t_quadraTrans_tilde *x = (t_quadraTrans_tilde *)(w[1]);
    t_sample  *inXL   =  (t_sample *)(w[2]); //inlet 1
    t_sample  *inXR   =  (t_sample *)(w[3]); //inlet 2
    t_sample  *outSR  =  (t_sample *)(w[4]); //outlet 4
    t_sample  *outSL  =  (t_sample *)(w[5]); //outlet 3
    t_sample  *outR   =  (t_sample *)(w[6]); //outlet 2
    t_sample  *outL   =  (t_sample *)(w[7]); //outlet 1
   
    int   blocksize   =         (int)(w[8]);

    int i;
    float mux = (float)1.0 / (float)x->fftsize;
    x->nbins = x->fftsize/2 + 1;

    float gain1=1;
    float gain2=0.5;
        
        for (i = 0; i < x->fftsize; i++) {
            if(i < blocksize){
                x->fftinL[i] = (inXL[i] + inXR[i])*gain1;
                x->fftinR[i] = (inXL[i] - inXR[i])*gain1;
            } else { // blocksize <= i <fftsize
                x->fftinL[i] = 0;
                x->fftinR[i] = 0;
            }
        }
        
        fftwf_execute(x->fftplanL);
        fftwf_execute(x->fftplanR);
        
        for (i = 0; i < x->fftsize; i++) {
            if(i < x->nPts){
                double tmp=(x->currentImpulse[0][i] + x->currentImpulse[1][i]);
                x->fftinHrirL[i] = 1 / (tmp==0?0.00000000000000000001:tmp);
                tmp=(x->currentImpulse[0][i] - x->currentImpulse[1][i]);
                x->fftinHrirR[i] = 1 / (tmp==0?0.00000000000000000001:tmp);
            } else { // nPts <= i <fftsize
                x->fftinHrirL[i] = 0;
                x->fftinHrirR[i] = 0;
            }
        }
        fftwf_execute(x->fftplanHrirL);
        fftwf_execute(x->fftplanHrirR);
        
        // convolve
        for (i = 0; i < x->nbins; i++) {
            x->fftinInvL[i][0] = (x->fftoutL[i][0] * x->fftoutHrirL[i][0] - x->fftoutL[i][1] * x->fftoutHrirL[i][1]);// * mux;
            x->fftinInvL[i][1] = (x->fftoutL[i][0] * x->fftoutHrirL[i][1] + x->fftoutL[i][1] * x->fftoutHrirL[i][0]);// * mux;
            x->fftinInvR[i][0] = (x->fftoutR[i][0] * x->fftoutHrirR[i][0] - x->fftoutR[i][1] * x->fftoutHrirR[i][1]);// * mux;
            x->fftinInvR[i][1] = (x->fftoutR[i][0] * x->fftoutHrirR[i][1] + x->fftoutR[i][1] * x->fftoutHrirR[i][0]);// * mux;
        }
        
        fftwf_execute(x->fftplanInvL);
        fftwf_execute(x->fftplanInvR);
        
        for(i = 0; i < x->fftsize; i++){

            float tmp;
            tmp = x->fftoutInvL[i]*mux;
            x->fftoutInvL[i] += x->fftoutInvR[i]*mux;
            x->fftoutInvR[i] = tmp - x->fftoutInvR[i]*mux;
            
            x->buffer[i][0] = x->buffer[i][0] + x->fftoutInvL[i];
            x->buffer[i][1] = x->buffer[i][1] + x->fftoutInvR[i];                                   

        }
        
        // Outputs
        int angle=(int)roundf(fabs(((float)x->aziL)))%360;
        float rad_ang=degToRad(angle);
    
        for (i = 0; i < x->fftsize; i++) {
            
            if(i < blocksize){   
                *outL++=((angle>315||angle<45?x->buffer[i][0]*cosf(rad_ang*2):0)+
                        (angle>270&&angle<360?(x->buffer[i][0]+x->buffer[i][1])*cosf(rad_ang*2-(2*315)):0)+
                        (angle>225&&angle<315?x->buffer[i][1]*cosf(rad_ang*2-(2*270)):0))*gain2;
                
                *outR++=((angle>315||angle<45?x->buffer[i][1]*cosf(rad_ang*2):0)+
                        (angle>0&&angle<90?(x->buffer[i][0]+x->buffer[i][1])*cosf(rad_ang*2-(2*45)):0)+
                        (angle>45&&angle<135?x->buffer[i][0]*cosf(rad_ang*2-(2*90)):0))*gain2;
                
                *outSL++=((angle>135&&angle<225?x->buffer[i][1]*cosf(rad_ang*2-(2*180)):0)+
                        (angle>180&&angle<270?(x->buffer[i][0]+x->buffer[i][1])*cosf(rad_ang*2-(2*225)):0)+
                        (angle>225&&angle<315?x->buffer[i][0]*cosf(rad_ang*2-(2*270)):0))*gain2;
                
                *outSR++=((angle>135&&angle<225?x->buffer[i][0]*cosf(rad_ang*2-(2*180)):0)+
                        (angle>90&&angle<180?(x->buffer[i][0]+x->buffer[i][1])*cosf(rad_ang*2-(2*135)):0)+
                        (angle>45&&angle<135?x->buffer[i][1]*cosf(rad_ang*2-(2*90)):0))*gain2;
            }
            x->buffer[i][0] = x->buffer[i+blocksize][0];
            x->buffer[i][1] = x->buffer[i+blocksize][1];
        }
    
    return (w+9);
}
    

void quadraTrans_tilde_dsp(t_quadraTrans_tilde *x, t_signal **sp){
    dsp_add(quadraTrans_tilde_perform, 8, x,
            sp[0]->s_vec, sp[1]->s_vec,
            sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
            sp[0]->s_n);
    
    int power, base;
    x->sr = sp[0]->s_sr;
    char file[2000] = "";
    char str[8] = "";
    
#ifdef WIN_VERSION
    strcpy_s(file, _countof(file), x->path);
    strcat_s(file, _countof(file), "/HRIR@");
    sprintf_s(str, sizeof(str), "%d", (int)x->sr);
    strcat_s(file, _countof(file), str);
    strcat_s(file, _countof(file), ".db");
    if (_access(file, F_OK) == -1) {
#endif
#ifdef PDMAC_VERSION
    strcpy(file,  x->path);
    strcat(file,  "/HRIR@");
    sprintf(str,  "%d", (int)x->sr);
    strcat(file,  str);
    strcat(file,  ".db");
    if (access(file, F_OK) == -1) {
#endif
        // file doesn't exist
        post("Warning: The database %s doesn't exist,\ntrying to use %s database", file, DBFILENAME);
            
#ifdef WIN_VERSION
        strcpy_s(file, _countof(file), x->path);
        strcat_s(file, _countof(file), "/HRIR@");
        sprintf_s(str, sizeof(file), "%d", 44100);
        strcat_s(file, _countof(file), str);
        strcat_s(file, _countof(file), ".db");
        if (_access(file, F_OK) == -1) {
#endif
#ifdef PDMAC_VERSION
        strcpy(file, x->path);
        strcat(file,  "/HRIR@");
        sprintf(str,  "%d", 44100);
        strcat(file,  str);
        strcat(file,  ".db");
        if (access(file, F_OK) == -1){
#endif
            post("Error: The database %s doesn't exist either, trying downloading from \n http://arts.u-aizu.ac.jp/research/hrir-with-distance-control/");
        } else {
            post("Warning: Using %s at a sampling freq. of %d, results are not correct.", file, (int)x->sr);
        }
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
        
        post("quadraTrans~: Max. blocksize: 8192, Max. taps %d, currently using %d \n", base, x->nPts);
    }
        
        x->convsize = x->nPts + sp[0]->s_n -1;
        x->fftsize = nextPo2(x->convsize);
        post("fftsize: %d, convsize: %d",x->fftsize,x->convsize);
        for(int i = 0; i < x->fftsize; i++){
            x->buffer[i][0] = 0.0;
            x->buffer[i][1] = 0.0;
        }
        
        //set L
        x->fftinL = fftwf_alloc_real(x->fftsize);
        x->fftoutL = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftplanL = fftwf_plan_dft_r2c_1d(x->fftsize, x->fftinL, x->fftoutL, FFTW_MEASURE);
        //set R
        x->fftinR = fftwf_alloc_real(x->fftsize);
        x->fftoutR = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftplanR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinR, x->fftoutR, FFTW_MEASURE);
        //set HrirL
        x->fftinHrirL = fftwf_alloc_real(x->fftsize);
        x->fftoutHrirL = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftplanHrirL = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirL, x->fftoutHrirL, FFTW_MEASURE);
        //set HrirR
        x->fftinHrirR = fftwf_alloc_real(x->fftsize);
        x->fftoutHrirR = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftplanHrirR = fftwf_plan_dft_r2c_1d((x->fftsize), x->fftinHrirR, x->fftoutHrirR, FFTW_MEASURE);
        //set invL
        x->fftinInvL = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftoutInvL = fftwf_alloc_real(x->fftsize);
        x->fftplanInvL = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvL, x->fftoutInvL, FFTW_MEASURE);
        //set invR
        x->fftinInvR = fftwf_alloc_complex(x->fftsize/2+1);
        x->fftoutInvR = fftwf_alloc_real(x->fftsize);
        x->fftplanInvR = fftwf_plan_dft_c2r_1d((x->fftsize), x->fftinInvR, x->fftoutInvR, FFTW_MEASURE);
        
        findFilter(x, 160, 0, 45);   // spearker L
        memcpy(&x->speackersFilters[0][0],x->currentImpulse[0],x->nPts * sizeof(float));
        memcpy(&x->speackersFilters[0][1],x->currentImpulse[0],x->nPts * sizeof(float));
        // findFilter(x, 160, 0, 135);  // speaker SL
        // memcpy(&x->speackersFilters[2][0],x->currentImpulse[0],x->nPts * sizeof(float));
        // memcpy(&x->speackersFilters[2][1],x->currentImpulse[0],x->nPts * sizeof(float));
        // findFilter(x, 160, 0, 225); // speaker SR
        // memcpy(&x->speackersFilters[3][0],x->currentImpulse[0],x->nPts * sizeof(float));
        // memcpy(&x->speackersFilters[3][1],x->currentImpulse[0],x->nPts * sizeof(float));
        // findFilter(x, 160, 0, 315);  // speaker R
        // memcpy(&x->speackersFilters[1][0],x->currentImpulse[0],x->nPts * sizeof(float));
        // memcpy(&x->speackersFilters[1][1],x->currentImpulse[0],x->nPts * sizeof(float));
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
    
    if (x->connected) {
        sqlite3_close(x->db);
        post("measurement database closed");
        x->connected = 0;
    }
    
    // fftw free
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
    
    //inlet
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->aziL);
    
    //outlet
    x->x_outL  = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outR  = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outSL = outlet_new(&x->x_obj, gensym("signal"));
    x->x_outSR = outlet_new(&x->x_obj, gensym("signal"));
    
    x->nPtsInDb = 300;
    x->nPts = DEFAULT_N_POINTS;
    
    strcat(x->path, canvas_getdir(canvas_getcurrent())->s_name);
    x->connected = 0;
    
    return (void *)x;
}

void quadraTrans_tilde_setup(void) {
    quadraTrans_tilde_class = class_new(gensym("quadraTrans~"),
                                        (t_newmethod)quadraTrans_tilde_new,
                                        (t_method)quadraTrans_tilde_free,
                                        sizeof(t_quadraTrans_tilde),
                                        CLASS_DEFAULT,
                                        A_DEFSYMBOL, A_DEFSYMBOL, A_DEFFLOAT, 0);
    
    class_addmethod(quadraTrans_tilde_class, (t_method)quadraTrans_tilde_dsp, gensym("dsp"), 0);
    CLASS_MAINSIGNALIN(quadraTrans_tilde_class, t_quadraTrans_tilde, f);
}
