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

#define round(x) ((x) < LONG_MIN-0.5 ||(x) > LONG_MAX+0.5 ?overflowError((x)):((x)>=0?(long)((x)+0.5):(long)((x)-0.5)))

/**
Used in round() for calculating the overflow error
*/
long overflowError(float y) {
    return (y>0.0) ? round(y - LONG_MAX) : round(y + LONG_MAX);
}

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

float Float(char unsigned *const p) {
    float val;
    memcpy(&val, p, sizeof val);
    return val;
}

/** This function computes the area A of a triangle using Heron's formula.
@param x, y, z, are the lengths of the triangle sides
*/
float triangleArea(float x, float y, float z) {
    return pow(((x + y + z)*(-x + y + z)*(x - y + z)*(x + y - z)) / 16, 0.5);
}

/** This function returns the euclidian distance between a given point q and a measurement in the DB.
@param x A pointer to the inner structure of the hrir~ object
@param i the current measurement being evaluated.
*/
float euclidianD(t_quadraTrans_tilde *x, int i) {
    float r, e, a;
    r = x->dis - x->set[i][0];
    e = x->ele - x->set[i][1];
    a = x->azi - x->set[i][2];
    return pow((r*r + e*e + a*a), 0.5);
}

/** This function returns the euclidian distance between 2 measurements in the DB.
@param x A pointer to the inner structure of the hrir~ object
@param i the first measurement
@param j the second measurement
*/
float euclidianD_2(t_quadraTrans_tilde *x, int i, int j) {
    float r, e, a;
    r = x->set[j][0] - x->set[i][0];
    e = x->set[j][1] - x->set[i][1];
    a = x->set[j][2] - x->set[i][2];
    return pow((r*r + e*e + a*a), 0.5);
}

/** This function computes the squared distance between the current point and four measurement points, returning the closest three in the "set" array of x.
This is a auxiliary function for the barimetric interpolation.
@param x A pointer to the inner structure of the hrir~ object
@param n The number of points considered for the interpolation
*/
int computeSquareDistance(t_quadraTrans_tilde *x, int n, int mixmax) {
    int i, minmaxID;
    float r, e, a, maxDist;
    minmaxID = (mixmax == 0) ? 99999999 : 0;
    maxDist = 0;
    for (i = 0; i < n; i++) {
        r = x->dis - x->set[i][0];
        e = x->ele - x->set[i][1];
        a = x->azi - x->set[i][2];
        x->dist[i] = (r*r + e*e + a*a);
        if (mixmax == 0) { // min
            if (x->dist[i] < maxDist) {
                maxDist = x->dist[i];
                minmaxID = i;
            }
        } else { // max
            if (x->dist[i] > maxDist) {
                maxDist = x->dist[i];
                minmaxID = i;
            }
        }

    }
    return minmaxID;
}

/** This function takes the most measurement and place it beyond the number of points needed for the interpolation
@param x A pointer to the inner structure of the hrir~ object
@param n The number of points considered for the interpolation
*/
void quicksort(t_quadraTrans_tilde *x, int n) {
    int i, maxId;
    maxId = computeSquareDistance(x, n, 1);
    x->set[n][0] = x->set[maxId][0];
    x->set[n][1] = x->set[maxId][1];
    x->set[n][2] = x->set[maxId][2];
    for (i = maxId; i < n; i++) {
        x->set[i][0] = x->set[i + 1][0];
        x->set[i][1] = x->set[i + 1][1];
        x->set[i][2] = x->set[i + 1][2];
    }
}

/** This function computes the volume V of a tetrahedron defined by points selected according to the index i
@param i An integer indicating the index to choose the points for the volume
*/
float theVolume(t_quadraTrans_tilde *x, int n) {
    float p[4][4];
    int i, j, k;
    float mult;
    float deter = 1;

    p[0][0] = x->dis;
    p[0][1] = x->ele;
    p[0][2] = x->azi;
    p[0][3] = 1;
    switch (n) {
    case 0: // [0 1 2]
        p[1][0] = x->set[0][0];
        p[1][1] = x->set[0][1];
        p[1][2] = x->set[0][2];
        p[1][3] = 1;
        p[2][0] = x->set[1][0];
        p[2][1] = x->set[1][1];
        p[2][2] = x->set[1][2];
        p[2][3] = 1;
        p[3][0] = x->set[2][0];
        p[3][1] = x->set[2][1];
        p[3][2] = x->set[2][2];
        p[3][3] = 1;
        break;
    case 1: // [0 1 3]
        p[1][0] = x->set[0][0];
        p[1][1] = x->set[0][1];
        p[1][2] = x->set[0][2];
        p[1][3] = 1;
        p[2][0] = x->set[1][0];
        p[2][1] = x->set[1][1];
        p[2][2] = x->set[1][2];
        p[2][3] = 1;
        p[3][0] = x->set[3][0];
        p[3][1] = x->set[3][1];
        p[3][2] = x->set[3][2];
        p[3][3] = 1;
        break;
    case 2: // [0 2 3]
        p[1][0] = x->set[0][0];
        p[1][1] = x->set[0][1];
        p[1][2] = x->set[0][2];
        p[1][3] = 1;
        p[2][0] = x->set[2][0];
        p[2][1] = x->set[2][1];
        p[2][2] = x->set[2][2];
        p[2][3] = 1;
        p[3][0] = x->set[3][0];
        p[3][1] = x->set[3][1];
        p[3][2] = x->set[3][2];
        p[3][3] = 1;
        break;
    case 3: // [1 2 3]
        p[1][0] = x->set[1][0];
        p[1][1] = x->set[1][1];
        p[1][2] = x->set[1][2];
        p[1][3] = 1;
        p[2][0] = x->set[2][0];
        p[2][1] = x->set[2][1];
        p[2][2] = x->set[2][2];
        p[2][3] = 1;
        p[3][0] = x->set[3][0];
        p[3][1] = x->set[3][1];
        p[3][2] = x->set[3][2];
        p[3][3] = 1;
        break;
    default:
        break;
    }

    for (i = 0; i<4; i++) {
        for (j = 0; j<4; j++) {
            mult = p[j][i] / p[i][i];
            for (k = 0; k<4; k++) {
                if (i == j)
                    break;
                p[j][k] = p[j][k] - p[i][k] * mult;
            }
        }
    }
    for (i = 0; i<4; i++) {
        deter = deter*p[i][i];
    }
    return fabs(deter) / 6;
}

static void makewindow(double *w, int n) {
    int i;
    double xshift = n / 2.0;
    double x;
    for (i = 0; i < n; i++) {
        x = (i - xshift) / xshift;
        w[i] = 0.5 * (1 + cos(MyPI * x));
    }
}

int nextPo2(t_float n){
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
