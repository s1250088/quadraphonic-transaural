/*
*  hrirFunctions.h
*
* This file comprises all the supplementary functions needed for the hrir~ object.
*
* See licence in the file hrir~.c
*
* Created by Julian Villegas on 2013.
* @author Julian Villegas <julian ^_^ at ^_^ u-aizu (*) ac (*) jp>
* @version 0.2a
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "m_pd.h"

//To compile it for Pd on Mac
# define PDMAC_VERSION


/**
calculates an integer value from a float value passed from the GUI (useful for azimuth)
*/
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

/** Validates the data passed to hrir~ inlets.
* Distances are between 20 and 160 cm,
* elevation between -40 and 90\[\degree\], and
* azimuths between 0 and 360\[\degree\].
* @param x A pointer to the inner structure of the hrir~ object.
*/
/*static void validatePoint(t_quadraTrans_tilde *x) {
    if (x->dis < 20) {
        x->dis = 20;
    }
    if (x->dis > 160) {
        x->dis = 160;
    }
    if (x->ele < -40)
        x->ele = -40;
    if (x->ele > 90)
        x->ele = 90;

    if (x->ele == 90) {
        x->azi = 0;
    } else {
        x->azi = (int)(x->azi) % 360;
        x->azi = (x->azi<0) ? x->azi + 360 : x->azi;
    }
}*/

/** finds the 2 closest measured distances from a given point
* @param x A pointer to the inner structure of the hrir~ object.
*/
static void rangeDis(t_quadraTrans_tilde *x) {
    t_float value;
    x->on[0] = -1;
    value = 10 * ((int)x->dis / 10);
    x->rs[0] = value;
    if (value < 160) {
        x->rs[1] = value + 10;
    } else {
        x->rs[1] = 160;
    }
    if ((x->rs[0] == x->dis) || (x->rs[1] == x->dis)) {
        x->on[0] = 1;
        x->rs[0] = x->dis;
    }
}

/** finds the 2 closest measured elevations from a given point
* @param x A pointer to the inner structure of the hrir~ object.
*/
static void rangeEle(t_quadraTrans_tilde *x) {
    t_float value;
    x->on[1] = -1;
    value = 10 * ((int)x->ele / 10);
    if (x->ele >= 0) {
        x->es[0] = value;
        if (value < 90) {
            x->es[1] = value + 10;
        } else {
            x->es[1] = 90;
        }
    } else {
        x->es[0] = value;
        if (value > -40) {
            x->es[1] = value - 10;
        } else {
            x->es[1] = -40;
        }
    }
    if ((x->es[0] == x->ele) || (x->es[1] == x->ele)) {
        x->on[1] = 1;
        x->es[0] = x->ele;
    }
}

/** finds the 2 closest measured azimuths from a given point
* @param x A pointer to the inner structure of the hrir~ object.
*/
static void rangeAzi(t_quadraTrans_tilde *x) {
    t_float value;
    x->on[2] = -1;
    value = 5 * ((int)x->azi / 5);
    x->as[0] = value;
    if (value < 360) {
        x->as[1] = value + 5;
    } else {
        x->as[1] = 360;
    }
    if ((x->as[0] == x->azi) || (x->as[1] == x->azi)) {
        x->on[2] = 1;
        x->as[0] = x->azi;
    }
}

/** form the candidate set of measurements for convolution
* @param x A pointer to the inner structure of the hrir~ object.
* to do: take care of the computation when the measurement is between 80 and 90 degrees of elevation
*/
static void formSet(t_quadraTrans_tilde *x) {
    int i, j, k;
    int counter = 0;
    int limR = 0, limE = 0, limA = 0;
    limR = (x->on[0]>0) ? 1 : 2;
    limE = (x->on[1]>0) ? 1 : 2;
    limA = (x->on[2]>0) ? 1 : 2;

    for (i = 0; i<limR; i++) {
        for (j = 0; j<limE; j++) {
            for (k = 0; k<limA; k++) {
                x->set[counter][0] = x->rs[i];
                x->set[counter][1] = x->es[j];
                x->set[counter][2] = x->as[k];
                counter += 1;
            }
        }
    }
    x->n_meas = counter;
}

/**
the db callback used to load into arrays the impulse response
*/
static int retrieveFromDd(t_quadraTrans_tilde *x, int argc, char **argv, char **azColName) {
    int i, j;
    float tmp[MAX_N_POINTS * 2];
    for (i = 0; i<argc; i++) {
        if (!strcmp(azColName[i], "hrir")) {
#ifdef WIN_VERSION
            memcpy_s(tmp, sizeof(tmp), argv[i], sizeof(float)*x->nPtsInDb);
#endif
#ifdef PDMAC_VERSION
            memcpy(tmp, argv[i], sizeof(float)*x->nPtsInDb);
#endif
            for (j = 0; j < x->nPts; j++) {
                x->currentImpulse[x->currentRow][0][j] = tmp[j];
                x->currentImpulse[x->currentRow][1][j] = tmp[j + x->nPtsInDb / 2];
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
#ifdef WIN_VERSION
    strcpy_s(query, _countof(query), "");
#endif
#ifdef PDMAC_VERSION
    strcpy(query, "");
#endif
#ifdef WIN_VERSION
    error = sprintf_s(query, sizeof(query), "SELECT hrir FROM measurements "
                      "WHERE r=%f and e=%f and a=%f "
                      , r, e, a);
#endif
#ifdef PDMAC_VERSION
    error = sprintf(query, "SELECT hrir FROM measurements "
                    "WHERE r=%f and e=%f and a=%f "
                    , r, e, a);
#endif
    if (error<0)
        post("error creating db query");
    else {
        if (x->connected) {
            x->rc = sqlite3_exec(x->db, query, retrieveFromDd, x, &(x->zErrMsg));
            checkError(x, "in retrieving db rows. ");
        }
    }
}

/** this  function computes the linear interpolation between a and b, given a
collinear point q. The interpolation is done via lerp. The points q, a, b are in (range,elevation,azimuth) form.
@param x A pointer to the inner structure of the hrir~ object
*/
static void linInterp(t_quadraTrans_tilde *x) {
    int i;
    t_float weight[3];
    //t_float itd[2];
    for (i = 0; i<2; i++) {
        x->currentRow = i;
        findFilter(x, x->set[i][0], x->set[i][1], x->set[i][2]);
        //itd[i] = x->currentItd;
        weight[i] = euclidianD(x, i);
    }
    weight[2] = weight[0] + weight[1];
    for (i = 0; i<x->nPts; i++) {
        x->currentImpulse[0][0][i] = weight[1] * x->currentImpulse[0][0][i];
        x->currentImpulse[0][1][i] = weight[1] * x->currentImpulse[0][1][i];
        x->currentImpulse[1][0][i] = weight[0] * x->currentImpulse[1][0][i];
        x->currentImpulse[1][1][i] = weight[0] * x->currentImpulse[1][1][i];
        x->currentImpulse[8][0][i] = (x->currentImpulse[0][0][i] + x->currentImpulse[1][0][i]) / weight[2];
        x->currentImpulse[8][1][i] = (x->currentImpulse[0][1][i] + x->currentImpulse[1][1][i]) / weight[2];
    }
    //x->currentItd = (weight[0]*itd[1] + weight[1]*itd[0]) / weight[2];
}

/** This function computes the barycentric interpolation between three points
given a coplanar point q. The interpolation is done via lerp.
The points are in (range,elevation,azimuth) form.
@param x A pointer to the inner structure of the hrir~ object
*/
static void baryInterp(t_quadraTrans_tilde *x) {
    int i;
    t_float distances[6]; ///distances from q to A, B, C, then AB, AC, BC.
    t_float weight[4]; ///ABQ, ACQ, BCQ, ABC.
    //t_float itd[3];
    quicksort(x, 4);

    for (i = 0; i < 3; i++) {
        x->currentRow = i;
        findFilter(x, x->set[i][0], x->set[i][1], x->set[i][2]);
        //itd[i] = x->currentItd;
        distances[i] = euclidianD(x, i);
    }
    distances[3] = euclidianD_2(x, 0, 1); //AB
    distances[4] = euclidianD_2(x, 0, 2); //AC
    distances[5] = euclidianD_2(x, 1, 2); //BC

    weight[0] = triangleArea(distances[0], distances[1], distances[3]);
    weight[1] = triangleArea(distances[0], distances[2], distances[4]);
    weight[2] = triangleArea(distances[1], distances[2], distances[5]);
    weight[3] = weight[0] + weight[1] + weight[2];

    for (i = 0; i<x->nPts; i++) {
        x->currentImpulse[0][0][i] = weight[2] * x->currentImpulse[0][0][i];
        x->currentImpulse[0][1][i] = weight[2] * x->currentImpulse[0][1][i];

        x->currentImpulse[1][0][i] = weight[1] * x->currentImpulse[1][0][i];
        x->currentImpulse[1][1][i] = weight[1] * x->currentImpulse[1][1][i];

        x->currentImpulse[2][0][i] = weight[0] * x->currentImpulse[2][0][i];
        x->currentImpulse[2][1][i] = weight[0] * x->currentImpulse[2][1][i];

        x->currentImpulse[8][0][i] = (x->currentImpulse[0][0][i] +
                                      x->currentImpulse[1][0][i] +
                                      x->currentImpulse[2][0][i]) / weight[3];
        x->currentImpulse[8][1][i] = (x->currentImpulse[0][1][i] +
                                      x->currentImpulse[1][1][i] +
                                      x->currentImpulse[2][1][i]) / weight[3];
    }

    //x->currentItd = (weight[0]*itd[2] + weight[1]*itd[1] + weight[2]*itd[0]) / weight[3];
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

/** This function computes the volumetric interpolation between 4 points
given a point q that is included. The interpolation is done via lerp.
The points are in (range,elevation,azimuth) form.
@param x A pointer to the inner structure of the hrir~ object
*/
static void simplexInterp(t_quadraTrans_tilde *x) {
    int i;
    t_float weight[5];
    //t_float itd[4];

    for (i = 0; i < 4; i++) {
        x->currentRow = i;
        findFilter(x, x->set[i][0], x->set[i][1], x->set[i][2]);
        //itd[i] = x->currentItd;
        weight[i] = theVolume(x, i);
    }
    weight[4] = weight[0] + weight[1] + weight[2] + weight[3];

    for (i = 0; i<x->nPts; i++) {
        x->currentImpulse[0][0][i] = weight[3] * x->currentImpulse[0][0][i];
        x->currentImpulse[0][1][i] = weight[3] * x->currentImpulse[0][1][i];

        x->currentImpulse[1][0][i] = weight[2] * x->currentImpulse[1][0][i];
        x->currentImpulse[1][1][i] = weight[2] * x->currentImpulse[1][1][i];

        x->currentImpulse[2][0][i] = weight[1] * x->currentImpulse[2][0][i];
        x->currentImpulse[2][1][i] = weight[1] * x->currentImpulse[2][1][i];

        x->currentImpulse[3][0][i] = weight[0] * x->currentImpulse[3][0][i];
        x->currentImpulse[3][1][i] = weight[0] * x->currentImpulse[3][1][i];

        x->currentImpulse[8][0][i] = (x->currentImpulse[0][0][i] +
                                      x->currentImpulse[1][0][i] +
                                      x->currentImpulse[2][0][i] +
                                      x->currentImpulse[3][0][i]) / weight[4];

        x->currentImpulse[8][1][i] = (x->currentImpulse[0][1][i] +
                                      x->currentImpulse[1][1][i] +
                                      x->currentImpulse[2][1][i] +
                                      x->currentImpulse[3][1][i]) / weight[4];
    }
    //x->currentItd = (weight[0]*itd[3] +
    //                 weight[1]*itd[2] +
    //                 weight[2]*itd[1] +
    //                 weight[3]*itd[0]) / weight[4];

}

/** This function select the closest 4 measurement around the point and invoke the interpolation between them. The points are indexed like:
% select simplex
%   6___7   e  r
%  /|  /|   | /
% 2_|_3 |   |/__a
% | | | |
% | 4_|_5
% |/  |/
% 0___1
@param x A pointer to the inner structure of the hrir~ object
*/
static void tetraInterp(t_quadraTrans_tilde *x) {
    int idx;
    idx = computeSquareDistance(x, 8, 0); //find the closest measurement
    switch (idx) {
    case 0: // use the [0 1 2 4] simplex
        x->set[3][0] = x->set[4][0];
        x->set[3][1] = x->set[4][1];
        x->set[3][2] = x->set[4][2];
        simplexInterp(x);
        break;
    case 1: // use the [1 0 3 5] simplex
        x->set[8][0] = x->set[0][0];
        x->set[8][1] = x->set[0][1];
        x->set[8][2] = x->set[0][2];

        x->set[0][0] = x->set[1][0];
        x->set[0][1] = x->set[1][1];
        x->set[0][2] = x->set[1][2];

        x->set[1][0] = x->set[8][0];
        x->set[1][1] = x->set[8][1];
        x->set[1][2] = x->set[8][2];

        x->set[2][0] = x->set[3][0];
        x->set[2][1] = x->set[3][1];
        x->set[2][2] = x->set[3][2];

        x->set[3][0] = x->set[5][0];
        x->set[3][1] = x->set[5][1];
        x->set[3][2] = x->set[5][2];
        simplexInterp(x);
        break;
    case 2: //use the [2 0 3 6] simplex
        x->set[1][0] = x->set[0][0];
        x->set[1][1] = x->set[0][1];
        x->set[1][2] = x->set[0][2];

        x->set[0][0] = x->set[2][0];
        x->set[0][1] = x->set[2][1];
        x->set[0][2] = x->set[2][2];

        x->set[2][0] = x->set[3][0];
        x->set[2][1] = x->set[3][1];
        x->set[2][2] = x->set[3][2];

        x->set[3][0] = x->set[6][0];
        x->set[3][1] = x->set[6][1];
        x->set[3][2] = x->set[6][2];
        simplexInterp(x);
        break;
    case 3: //use the [3 1 2 7] simplex
        x->set[0][0] = x->set[3][0];
        x->set[0][1] = x->set[3][1];
        x->set[0][2] = x->set[3][2];

        x->set[3][0] = x->set[7][0];
        x->set[3][1] = x->set[7][1];
        x->set[3][2] = x->set[7][2];
        simplexInterp(x);
        break;
    case 4: //use [4 0 5 6]
        x->set[1][0] = x->set[0][0];
        x->set[1][1] = x->set[0][1];
        x->set[1][2] = x->set[0][2];

        x->set[0][0] = x->set[4][0];
        x->set[0][1] = x->set[4][1];
        x->set[0][2] = x->set[4][2];

        x->set[2][0] = x->set[5][0];
        x->set[2][1] = x->set[5][1];
        x->set[2][2] = x->set[5][2];

        x->set[3][0] = x->set[6][0];
        x->set[3][1] = x->set[6][1];
        x->set[3][2] = x->set[6][2];
        simplexInterp(x);
        break;
    case 5: // use [5 1 4 7]
        x->set[0][0] = x->set[5][0];
        x->set[0][1] = x->set[5][1];
        x->set[0][2] = x->set[5][2];

        x->set[2][0] = x->set[4][0];
        x->set[2][1] = x->set[4][1];
        x->set[2][2] = x->set[4][2];

        x->set[3][0] = x->set[7][0];
        x->set[3][1] = x->set[7][1];
        x->set[3][2] = x->set[7][2];
        simplexInterp(x);
        break;
    case 6: // use [6 2 4 7]
        x->set[0][0] = x->set[6][0];
        x->set[0][1] = x->set[6][1];
        x->set[0][2] = x->set[6][2];

        x->set[1][0] = x->set[2][0];
        x->set[1][1] = x->set[2][1];
        x->set[1][2] = x->set[2][2];

        x->set[2][0] = x->set[4][0];
        x->set[2][1] = x->set[4][1];
        x->set[2][2] = x->set[4][2];

        x->set[3][0] = x->set[7][0];
        x->set[3][1] = x->set[7][1];
        x->set[3][2] = x->set[7][2];
        simplexInterp(x);
        break;
    case 7: // use [7 3 5 6]
        x->set[0][0] = x->set[7][0];
        x->set[0][1] = x->set[7][1];
        x->set[0][2] = x->set[7][2];

        x->set[1][0] = x->set[3][0];
        x->set[1][1] = x->set[3][1];
        x->set[1][2] = x->set[3][2];

        x->set[2][0] = x->set[5][0];
        x->set[2][1] = x->set[5][1];
        x->set[2][2] = x->set[5][2];

        x->set[3][0] = x->set[6][0];
        x->set[3][1] = x->set[6][1];
        x->set[3][2] = x->set[6][2];
        simplexInterp(x);
        break;
    default:
        printf("Invalid\n");
    }

}


/** This function calls the adequated interpolation depending on the independence of coordinates between the point in question and the measurements.
* @param x A pointer to the inner structure of the hrir~ object.
*/
static void testCoplan_lin(t_quadraTrans_tilde *x) {
    switch (x->n_meas) {
    case 1:
        x->currentRow = 8;
        findFilter(x, x->set[0][0], x->set[0][1], x->set[0][2]);
        break;
    case 2:
        linInterp(x);
        break;
    case 4:
        baryInterp(x);
        break;
    case 8:
        tetraInterp(x);
        break;
    default:
        break;
    }
}


/**
Calculates the convolution between two arrays (the signal and the hrir)
*/
/*static void convolve(t_quadraTrans_tilde *x, int blocksize, t_float *audio_in, t_float *audio_outL, t_float *audio_outR) {
    float convSumI, convSumC; /// convolution sum of the ipsilateral and contralateral sides
    float inSample;
    int blockScale = MAX_BLOCKSIZE / blocksize;
    unsigned scaledBlocksize;
    unsigned blocksizeDelta;
    int i;

    while (blocksize--) {
        convSumI = 0;
        convSumC = 0;
        inSample = *(audio_in++);
        x->convBufferC[x->bufferPin] = inSample; // Contralateral ear
        x->convBufferI[x->bufferPin] = inSample; // Ipsilateral ear
        scaledBlocksize = blocksize * blockScale;
        blocksizeDelta = MAX_BLOCKSIZE - 1 - scaledBlocksize;

        for (i = 0; i < x->nPts; i++) {
            convSumC += x->convBufferC[(x->bufferPin - i) & (x->nPts - 1)] *
                        (x->currentImpulse[8][0][i] * x->crossCoef[scaledBlocksize] +
                         x->previousImpulse[0][i] * x->crossCoef[blocksizeDelta]);
            convSumI += x->convBufferI[(x->bufferPin - i) & (x->nPts - 1)] *
                        (x->currentImpulse[8][1][i] * x->crossCoef[scaledBlocksize] +
                         x->previousImpulse[1][i] * x->crossCoef[blocksizeDelta]);
            x->previousImpulse[0][i] = x->currentImpulse[8][0][i];
            x->previousImpulse[1][i] = x->currentImpulse[8][1][i];
        }
        x->bufferPin = (x->bufferPin + 1) & (x->nPts - 1);
        if (x->azi <= 180) { // The db only has the right hemisphere hrir
            *audio_outL++ = convSumC / 3; // 3 is to make it less than 1 at 20cm, 90 deg, 0 elevation when a white noise is being directionalized
            *audio_outR++ = convSumI / 3;
        } else {
            *audio_outR++ = convSumC / 3;
            *audio_outL++ = convSumI / 3;
        }
    }
}*/
