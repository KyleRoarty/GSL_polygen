#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmodule.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

#include "gslsg.h"

//Vertex dimensions. Generally 3, could be 2
#define V_DIM 3

#define DELTA 0.0000000001
#define SIGMA 0.00000000000000000000000000001

void print_vert_3(gsl_vector **vert, int num_v){
    for(int i = 0; i < num_v; i++){
        printf("Vertex %d\n", i);
        gsl_vector_fprintf(stdout, vert[i], "%f");
        printf("\n");
    }
}

void print_seg_3(seg_3 **seg, int idx){
    printf("Segment %d\nStart:\n", idx);
    gsl_vector_fprintf(stdout, seg[idx]->start, "%f");
    printf("End:\n");
    gsl_vector_fprintf(stdout, seg[idx]->end, "%f");
    printf("\n");
}

void print_seg_3_all(seg_3 **seg, int num_s){
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

//[from, to)
int sum(int from, int to){
    int ret = 0;
    for(int i = from; i > to; i--)
        ret += i;
    return ret;
}

// x < y; sz = num points (cla)
// x, y are points that define a segment
int segFromI(int x, int y, int sz){
    return x*sz-sum(x, 0)+y-(x+1);
}

gboolean print_g_tree(gpointer key, gpointer value, gpointer data){
    if(value == NULL)
        return true;

    printf("%2d, ", GPOINTER_TO_INT(value));
    return false;
}

int sort_integers(gconstpointer a, gconstpointer b){
    int int_a, int_b;
    int_a = GPOINTER_TO_INT(a);
    int_b = GPOINTER_TO_INT(b);

    return int_a-int_b;
}

//For a *_3 function, constants, like for gsl_vector_malloc, will be 3
int segment_slope_3(seg_3 *seg){
    int err;
    seg->slope = gsl_vector_calloc(3);

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

    free(num);
    return GSL_SUCCESS;
}

// TODO: replace V_DIM with a different constant here
int segment_intersect_3(seg_3 *seg1, seg_3 *seg2, gsl_vector *x){
    int err, ret;
    gsl_matrix *A;
    gsl_vector *tau, *b, *residual;

    residual = gsl_vector_alloc(V_DIM);
    tau = gsl_vector_alloc(V_DIM-1);
    b = gsl_vector_alloc(V_DIM);
    A = gsl_matrix_alloc(V_DIM, V_DIM-1);

    for(int i = 0; i < V_DIM; i++)
        gsl_vector_set(b, i, gsl_vector_get(seg2->start, i)-gsl_vector_get(seg1->start, i));

    for(int i = 0; i < V_DIM; i++){
        gsl_matrix_set(A, i, 0, gsl_vector_get(seg1->slope, i));
        gsl_matrix_set(A, i, 1, -gsl_vector_get(seg2->slope, i));
    }

    if( (err = gsl_linalg_QR_decomp(A, tau)) ){
        printf("Error: QR decomp. %d\n", err);
        ret = GSL_FAILURE;
        goto cleanup;
    }

    if( (err = gsl_linalg_QR_lssolve(A, tau, b, x, residual)) ){
        printf("Error: QR lssolve. %d\n", err);
        ret = GSL_FAILURE;
        goto cleanup;
    }

    gsl_vector_mul(residual, residual);
    if( gsl_vector_max(residual) > SIGMA ){
        ret = GSL_FAILURE;
        goto cleanup;
    }

    // TODO: if solution, check bounds for seg1, seg2 begin and end

    ret = GSL_SUCCESS;

cleanup:
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(tau);
    gsl_vector_free(residual);
    return ret;
}

// TODO: Double for loop; don't repeat code for seg1, seg2
int overlap_in_bounds_3(seg_3 *seg1, seg_3 *seg2, gsl_vector *x){

    gsl_vector *sum;
    int ret;

    sum = gsl_vector_calloc(seg1->start->size);

    //  Check if x is infinity, return -1 if it is
    for(int i = 0 ; i < x->size; i++){
        if(gsl_isinf(gsl_vector_get(x, i))){
            ret = GSL_FAILURE;
            goto cleanup;
        }
    }


    // slope*x+start
    gsl_vector_memcpy(sum, seg1->slope);
    gsl_vector_scale(sum, gsl_vector_get(x, 0));
    gsl_vector_add(sum, seg1->start);

    //gsl_vector_fprintf(stdout, sum, "%1.50f");

    //Do I need to do all of the checks?
    //Do the check if the start/end aren't equal
    for(int i = 0; i < 3; i++){
        if( gsl_vector_get(seg1->start, i) != gsl_vector_get(seg1->end, i) ){
            if( (gsl_vector_get(sum, i) < gsl_vector_get(seg1->start, i) && gsl_vector_get(sum, i) < gsl_vector_get(seg1->end, i)) ||
                (gsl_vector_get(sum, i) > gsl_vector_get(seg1->start, i) && gsl_vector_get(sum, i) > gsl_vector_get(seg1->end, i)) ){
                ret = GSL_FAILURE;
                goto cleanup;
            }
            // Close enough to the actual value
            if( fabs(gsl_vector_get(sum, i)-gsl_vector_get(seg1->start, i)) < DELTA ||
                fabs(gsl_vector_get(sum, i)-gsl_vector_get(seg1->end, i)) < DELTA ){
                ret = GSL_FAILURE;
                goto cleanup;
            }
        }
    }

    //gsl_vector_fprintf(stdout, sum, "%1.50f");

    gsl_vector_memcpy(sum, seg2->slope);
    gsl_vector_scale(sum, gsl_vector_get(x, 1));
    gsl_vector_add(sum, seg2->start);

    for(int i = 0; i < 3; i++){
        if( !(gsl_vector_get(seg2->start, i) == gsl_vector_get(seg2->end, i)) ){
            if( (gsl_vector_get(sum, i) < gsl_vector_get(seg2->start, i) && gsl_vector_get(sum, i) < gsl_vector_get(seg2->end, i)) ||
                (gsl_vector_get(sum, i) > gsl_vector_get(seg2->start, i) && gsl_vector_get(sum, i) > gsl_vector_get(seg2->end, i)) ){
                ret = GSL_FAILURE;
                goto cleanup;
            }
            // Close enough to the actual value
            if( fabs(gsl_vector_get(sum, i)-gsl_vector_get(seg2->start, i)) < DELTA ||
                fabs(gsl_vector_get(sum, i)-gsl_vector_get(seg2->end, i)) < DELTA ){
                ret = GSL_FAILURE;
                goto cleanup;
            }
        }
    }

    ret = GSL_SUCCESS;

cleanup:
    gsl_vector_free(sum);
    return ret;
}

// TODO: Don't repeat code for a, b;
// TODO: replace conditionals with case statement?
int resolve_overlap(GTree *safe, GTree *ignore, int a, int b){
    // For a, b: Check if in safe, in overlap
    bool safe_a=false, ignore_a=false, safe_b=false, ignore_b=false;

    if(g_tree_lookup(safe, GINT_TO_POINTER(a)))
        safe_a = true;
    if(g_tree_lookup(ignore, GINT_TO_POINTER(a)))
        ignore_a = true;

    if(g_tree_lookup(safe, GINT_TO_POINTER(b)))
        safe_b = true;
    if(g_tree_lookup(ignore, GINT_TO_POINTER(b)))
        ignore_b = true;

    //printf("s_a: %d, i_a: %d, s_b: %d, i_b: %d\n", safe_a, ignore_a, safe_b, ignore_b);

    if((safe_a && ignore_a) || (safe_b && ignore_b))
        return -1;

    if((safe_a && ignore_b) || (ignore_a && safe_b) || (ignore_a && ignore_b))
        return 0;

    if(safe_a && safe_b){
        g_tree_remove(safe, GINT_TO_POINTER(b));
        g_tree_insert(ignore, GINT_TO_POINTER(b), GINT_TO_POINTER(b));
        return 0;
    }

    if(safe_a){
        g_tree_insert(ignore, GINT_TO_POINTER(b), GINT_TO_POINTER(b));
        return 0;
    }

    if(safe_b){
        g_tree_insert(ignore, GINT_TO_POINTER(a), GINT_TO_POINTER(a));
        return 0;
    }

    if(ignore_a){
        g_tree_insert(safe, GINT_TO_POINTER(b), GINT_TO_POINTER(b));
        return 0;
    }

    if(ignore_b){
        g_tree_insert(safe, GINT_TO_POINTER(a), GINT_TO_POINTER(a));
        return 0;
    }

    g_tree_insert(ignore, GINT_TO_POINTER(a), GINT_TO_POINTER(a));
    g_tree_insert(safe, GINT_TO_POINTER(b), GINT_TO_POINTER(b));
    return 0;
}

gboolean same_point(gpointer key, gpointer value, gpointer data){
    args *comp_args;
    seg_3 *seg;
    comp_args = (args *)data;
    seg = comp_args->segs[GPOINTER_TO_INT(value)];
    if(seg->start == comp_args->comp || seg->end == comp_args->comp)
        *(comp_args->count) = *(comp_args->count) + 1;
    return false;
}

int main(int argc, char **argv){
    FILE *fp;
    //number of points/vertices, number of segments
    int num_v, num_s, num_t;

    //used to initialize segments
    int iter;

    //Check error from functions
    int err;

    //Used when counting overlapping overlaps (because that makes sense)
    int count;

    // TODO: Rename tri_ol to something that makes more sense
    args *tri_ol;

    gsl_vector **vert;
    gsl_vector *x;
    seg_3 **seg;
    GTree *safe, *ignore;

    int currTri;
    tri_3 **tri;

    if(argc != 3){
        //n is number of vertices
        //file is csv of 3d vertices
        printf("USAGE: shapegen <n> <file>\n");
        exit(-1);
    }

    tri_ol = malloc(sizeof(args));

    safe = g_tree_new(sort_integers);
    ignore = g_tree_new(sort_integers);

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

    for(int i = 0; i < num_s; i++){
        for(int j = i+1; j < num_s; j++){
            x = gsl_vector_calloc(V_DIM-1);
            if( segment_intersect_3(seg[i], seg[j], x) == GSL_SUCCESS ){
                if( overlap_in_bounds_3(seg[i], seg[j], x) == GSL_SUCCESS ){
                    printf("Seg %d and seg %d\n", i, j);
                    resolve_overlap(safe, ignore, i, j);
                }
            }
            gsl_vector_free(x);
        }
    }

    printf("Safe: \n");
    g_tree_foreach(safe, print_g_tree, GINT_TO_POINTER(52));
    printf("\n");
    printf("Ignore: \n");
    g_tree_foreach(ignore, print_g_tree, GINT_TO_POINTER(52));
    printf("\n");


    //Each triangle needs 3 unique segments
    //However the entire search space isn't covered
    num_t = (num_v)*(num_v-1)*(num_v-2)/6;
    printf("num_t: %d\n", num_t);
    num_t -= g_tree_nnodes(ignore)*(num_t*3/num_s);
    printf("num_t: %d\n", num_t);

    tri_ol->segs = seg;
    for(int i = 0; i < num_v; i++){
        count = 0;
        tri_ol->count = &count;
        tri_ol->comp = vert[i];
        g_tree_foreach(ignore, same_point, (gpointer)tri_ol);
        num_t += ((count)*(count-1)/2);
        printf("Count %d: %d\n", i, count);
    }

    printf("Number of triangles: %d\n", num_t);

    tri = malloc(num_t*sizeof(tri_3 *));
    currTri = 0;
    // i < num_v-2 maybe?
    for(int i = 0; i < num_v; i++){
        for(int j = i+1; j < num_v; j++){
            //Check if segment is in the no-go list
            if(g_tree_lookup(ignore, GINT_TO_POINTER(segFromI(i, j, num_v))))
                continue;
            for(int k = j+1; k < num_v; k++){
                //Check if i-k, j-k is in the no-go list
                if(g_tree_lookup(ignore, GINT_TO_POINTER(segFromI(i, k, num_v))) ||
                   g_tree_lookup(ignore, GINT_TO_POINTER(segFromI(j, k, num_v))))
                    continue;

                tri[currTri] = malloc(sizeof(tri_3));
                tri[currTri]->vert[0] = vert[i];
                tri[currTri]->idx[0] = i;
                tri[currTri]->vert[1] = vert[j];
                tri[currTri]->idx[1] = j;
                tri[currTri]->vert[2] = vert[k];
                tri[currTri]->idx[2] = k;
                currTri++;
            }
        }
    }

    for(int i = 0; i < num_t; i++)
        printf("Tri %d: (%d, %d, %d)\n", i, tri[i]->idx[0], tri[i]->idx[1], tri[i]->idx[2]);

    //Memory management section

    for(int i = num_s-1; i >= 0; i--){
        gsl_vector_free(seg[i]->slope);
        free(seg[i]);
    }

    for(int i = num_v-1; i >= 0; i--)
        gsl_vector_free(vert[i]);

    free(seg);
    free(vert);
    fclose(fp);

    g_tree_destroy(ignore);
    g_tree_destroy(safe);

    free(tri_ol);
    return(0);
}
