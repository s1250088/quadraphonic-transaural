#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "m_pd.h"

//To compile it for Pd on Mac
#define PDMAC_VERSION

#ifdef WIN_VERSION
#include <io.h>
#endif

/**
This function returns 0 if the argument is not a power of 2
*/
int isPowerOfTwo(unsigned int x) {
    return ((x != 0) && ((x & (~x + 1)) == x));
}

/**
Check and display errors, freeing the memory pointed by the sqlite3 error pointer
*/
void checkError(t_quadraTrans_tilde *x, char *msg) {
    if (x->rc != SQLITE_OK) {
        post("%s. SQL error: %s\n", msg, x->zErrMsg);
        sqlite3_free(x->zErrMsg);
    }
}

/**
 the db callback used to load into arrays the impulse response
 */
static int retrieveFromDd(t_quadraTrans_tilde *x, int argc, char **argv, char **azColName) {
    int i,j;
    float tmp[x->nPtsInDb * sizeof(float)];
    for(i=0; i<argc; i++) {
        if (!strcmp(azColName[i],"hrir")) {
            memcpy(&tmp, argv[i], x->nPtsInDb * sizeof(float));
            for (j=0; j < x->nPts; j++) {
                x->currentImpulse[0][j]=tmp[j];
                x->currentImpulse[1][j]=tmp[j + x->nPtsInDb/2];
            }
        }
    }
    return 0;
}

/**
 Find a measurement in the database since only the right ear measurements are stored, angles greater than 180 are transformed.
 * r,e,a, must be valid measures
 */
static void findFilter(t_quadraTrans_tilde *x, float r, float e, float a) {
    int error;
    char query[2000];
    if (a > 180) {
        a = 360 - a;
    }
    if (a == 360) {
        a = 0;
    }
    if (e == 90) {
        a = 0;
    }

    strcpy(query,"");
    error=sprintf(query,"SELECT hrir FROM measurements "
                  "WHERE r=%f and e=%f and a=%f "
                  , r, e, a);
    if(error<0)
        post("error creating db query");
    else {
        if(x->connected) {
            x->rc = sqlite3_exec(x->db, query, retrieveFromDd, x, &(x->zErrMsg));
            checkError(x,"in retrieving db rows. ");
        }
    }
}

int nextPo2(int n){
    int i=0, result;
    
    while(1){
        if(n <= powf(2.0, (float)i)){
            result = (int)powf(2, i);
            break;
        }
        i++;
    }
    return result;
}
