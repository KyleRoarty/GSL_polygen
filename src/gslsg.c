#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include "gslsg.h"

//Vertex dimensions. Generally 3, could be 2
#define V_DIM 3

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
    //number of points/vertices
    int num_v;
    gsl_vector **vert;

    if(argc != 3){
        //n is number of vertices
        //file is csv of 3d vertices
        printf("USAGE: shapegen <n> <file>\n");
        exit(-1);
    }

    num_v = atoi(argv[1]);
    vert = malloc(num_v*sizeof(gsl_vector *));

    fp = fopen(argv[2], "r");

    for(int i = 0; i < num_v; i++){
        vert[i] = gsl_vector_alloc(V_DIM);
        if(V_DIM == 3)
            if(create_vertex_3(fp, vert[i]))
                exit(-1);
    }

    for(int i = 0; i < num_v; i++){
        printf("Vertex %d\n", i);
        gsl_vector_fprintf(stdout, vert[i], "%f");
        printf("\n");
    }

    for(int i = num_v-1; i >= 0; i--)
        gsl_vector_free(vert[i]);


    free(vert);
    return(0);
}
