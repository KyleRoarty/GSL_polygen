#ifndef __GSLSG_H__
#define __GSLSG_H__

#include <gsl/gsl_vector.h>

typedef struct{
    gsl_vector *start;
    gsl_vector *end;
    gsl_vector *slope;
} seg_3;

typedef struct{
    gsl_vector *vert[3];
    int idx[3];
} tri_3;

typedef struct{
    seg_3 **segs;
    gsl_vector *comp;
    int *count;
} args;



#endif
