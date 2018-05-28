#ifndef __GSLSG_H__
#define __GSLSG_H__

#include <gmodule.h>
#include <gsl/gsl_vector.h>
#include <stdbool.h>

typedef struct{
    gsl_vector *vert[2];
    gsl_vector *slope;
    int idx[2];
} seg_3;

typedef struct{
    gsl_vector *vert[3];
    int idx[3];
    bool use;
    bool ignore;
} tri_3;

typedef struct{
    int n;
    int *tri;
} tri_share;

typedef struct{
    seg_3 **segs;
    gsl_vector *comp;
    int *count;
} args;

typedef struct{
    seg_3 **segs;
    GTree *tree;
    bool *at_key;
    int *lc;
    int *gc;
    int num_v;
    gpointer key;
} d_args;

#endif
