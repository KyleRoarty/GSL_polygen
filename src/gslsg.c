#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include "gslsg.h"

//Vertex dimensions. Generally 3, could be 2
#define V_DIM 3

void print_vert_3(gsl_vector **vert, int num_v){
    for(int i = 0; i < num_v; i++){
        printf("Vertex %d\n", i);
        gsl_vector_fprintf(stdout, vert[i], "%f");
        printf("\n");
    }
}

void print_seg_3(seg_3 **seg, int num_s){
    for(int i = 0; i < num_s; i++){
        printf("Segment %d\nStart:\n", i);
        gsl_vector_fprintf(stdout, seg[i]->start, "%f");
        printf("End:\n");
        gsl_vector_fprintf(stdout, seg[i]->end, "%f");
        printf("Slope:\n");
        gsl_vector_fprintf(stdout, seg[i]->slope, "%f");
        printf("\n");
    }
}

//For a *_3 function, constants, like for gsl_vector_malloc, will be 3
int segment_slope_3(seg_3 *seg){
    int err;
    seg->slope = gsl_vector_alloc(3);
    if( (err = gsl_vector_add(seg->slope, seg->end)) ){
        // Handle error
        return GSL_FAILURE;
    }
    if( (err = gsl_vector_sub(seg->slope, seg->start)) ){
        // Handle error
        return GSL_FAILURE;
    }

    return GSL_SUCCESS;
}

int create_vertex_3(FILE *fp, gsl_vector *v){
    char *num=NULL, *saveptr=NULL;
    size_t len;

    if(getline(&num, &len, fp) == -1){
        return GSL_FAILURE;
    }

    gsl_vector_set(v, 0, atof(strtok_r(num, ",", &saveptr)));
    gsl_vector_set(v, 1, atof(strtok_r(NULL, ",", &saveptr)));
    gsl_vector_set(v, 2, atof(saveptr));

    return GSL_SUCCESS;
}

int main(int argc, char **argv){
    FILE *fp;
    //number of points/vertices, number of segments
    int num_v, num_s;

    //used to initialize segments
    int iter;

    //Check error from functions
    int err;

    gsl_vector **vert;
    seg_3 **seg;

    if(argc != 3){
        //n is number of vertices
        //file is csv of 3d vertices
        printf("USAGE: shapegen <n> <file>\n");
        exit(-1);
    }

    num_v = atoi(argv[1]);
    vert = malloc(num_v*sizeof(gsl_vector *));

    //Num lines: nCr(num_v, 2) = num_v!/(2!*(num-2)!), num*(num-1)/2
    num_s = num_v*(num_v-1)/2;
    seg = malloc(num_s*sizeof(seg_3 *));

    fp = fopen(argv[2], "r");

    for(int i = 0; i < num_v; i++){
        vert[i] = gsl_vector_alloc(V_DIM);
        if(V_DIM == 3)
            if(create_vertex_3(fp, vert[i]))
                exit(-1);
    }

    print_vert_3(vert, num_v);

    iter = 0;
    for(int i = 0; i < num_v; i++){
        for(int j = i+1; j < num_v; j++){
            seg[iter] = malloc(sizeof(seg_3));
            seg[iter]->start = vert[i];
            seg[iter]->end = vert[j];
            if( (err = segment_slope_3(seg[iter])) ){
                //handle error
                printf("Error! From segment_slope_3\n");
            }
            iter++;
        }
    }

    print_seg_3(seg, num_s);

    //Memory management section

    for(int i = num_s-1; i >= 0; i--){
        gsl_vector_free(seg[i]->slope);
        free(seg[i]);
    }

    for(int i = num_v-1; i >= 0; i--)
        gsl_vector_free(vert[i]);

    free(seg);
    free(vert);
    return(0);
}
