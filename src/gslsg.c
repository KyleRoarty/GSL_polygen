//TODO: Error when things break, don't just keep progressing

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
    gsl_vector_fprintf(stdout, seg[idx]->vert[0], "%f");
    printf("End:\n");
    gsl_vector_fprintf(stdout, seg[idx]->vert[1], "%f");
    printf("\n");
}

void print_seg_3_all(seg_3 **seg, int num_s){
    for(int i = 0; i < num_s; i++){
        printf("Segment %d\nStart:\n", i);
        gsl_vector_fprintf(stdout, seg[i]->vert[0], "%f");
        printf("End:\n");
        gsl_vector_fprintf(stdout, seg[i]->vert[1], "%f");
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

    seg_3 **seg = (seg_3 **)data;
    int i = GPOINTER_TO_INT(value)-1;

    printf("%2d (%d, %d); ", i, seg[i]->idx[0], seg[i]->idx[1]);
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

    if( (err = gsl_vector_add(seg->slope, seg->vert[1])) ){
        // Handle error
        return GSL_FAILURE;
    }
    if( (err = gsl_vector_sub(seg->slope, seg->vert[0])) ){
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
        gsl_vector_set(b, i, gsl_vector_get(seg2->vert[0], i)-gsl_vector_get(seg1->vert[0], i));

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

    ret = GSL_SUCCESS;

cleanup:
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(tau);
    gsl_vector_free(residual);
    return ret;
}

int overlap_in_bounds_3(seg_3 *seg1, seg_3 *seg2, gsl_vector *x){

    seg_3 *segs[2];
    gsl_vector *sum;
    int ret;

    segs[0] = seg1;
    segs[1] = seg2;

    sum = gsl_vector_calloc(segs[0]->vert[0]->size);

    //  Check if x is infinity, return -1 if it is
    for(int i = 0 ; i < x->size; i++){
        if(gsl_isinf(gsl_vector_get(x, i))){
            ret = GSL_FAILURE;
            goto cleanup;
        }
    }

    for(int j = 0; j < 2; j++){
        // slope*x+start
        gsl_vector_memcpy(sum, segs[j]->slope);
        gsl_vector_scale(sum, gsl_vector_get(x, 0));
        gsl_vector_add(sum, segs[j]->vert[0]);

        //gsl_vector_fprintf(stdout, sum, "%1.50f");

        //Do I need to do all of the checks?
        //Do the check if the start/end aren't equal
        for(int i = 0; i < 3; i++){
            if( gsl_vector_get(segs[j]->vert[0], i) != gsl_vector_get(segs[j]->vert[1], i) ){
                if( (gsl_vector_get(sum, i) < gsl_vector_get(segs[j]->vert[0], i) && gsl_vector_get(sum, i) < gsl_vector_get(segs[j]->vert[1], i)) ||
                    (gsl_vector_get(sum, i) > gsl_vector_get(segs[j]->vert[0], i) && gsl_vector_get(sum, i) > gsl_vector_get(segs[j]->vert[1], i)) ){
                    ret = GSL_FAILURE;
                    goto cleanup;
                }
                // Close enough to the actual value
                if( fabs(gsl_vector_get(sum, i)-gsl_vector_get(segs[j]->vert[0], i)) < DELTA ||
                    fabs(gsl_vector_get(sum, i)-gsl_vector_get(segs[j]->vert[1], i)) < DELTA ){
                    ret = GSL_FAILURE;
                    goto cleanup;
                }
            }
        }
    }
    ret = GSL_SUCCESS;

cleanup:
    gsl_vector_free(sum);
    return ret;
}

int resolve_overlap(GTree *safe, GTree *ignore, int a, int b){
    // For a, b: Check if in safe, in overlap
    int idx[2] = {a, b};
    bool safe_i[2] = {false, false};
    bool ignore_i[2] = {false, false};

    for(int i = 0; i < 2; i++){
        if(g_tree_lookup(safe, GINT_TO_POINTER(idx[i])))
            safe_i[i] = true;
        if(g_tree_lookup(ignore, GINT_TO_POINTER(idx[i])))
            ignore_i[i] = true;
    }

    //printf("s_a: %d, i_a: %d, s_b: %d, i_b: %d\n", safe_i[0], ignore_i[0], safe_i[1], ignore_i[1]);

    if((safe_i[0] && ignore_i[0]) || (safe_i[1] && ignore_i[1]))
        return -1;

    if((safe_i[0] && ignore_i[1]) || (ignore_i[0] && safe_i[1]) || (ignore_i[0] && ignore_i[1]))
        return 0;

    if(safe_i[0] && safe_i[1]){
        g_tree_remove(safe, GINT_TO_POINTER(idx[1]));
        g_tree_insert(ignore, GINT_TO_POINTER(idx[1]), GINT_TO_POINTER(idx[1]+1));
        return 0;
    }

    for(int i = 0; i < 2; i++){
        if(safe_i[i]){
            g_tree_insert(ignore, GINT_TO_POINTER(idx[i]), GINT_TO_POINTER(idx[i]+1));
            return 0;
        }
        if(ignore_i[i]){
            g_tree_insert(safe, GINT_TO_POINTER(idx[i+1 % 2]), GINT_TO_POINTER(idx[i+1 % 2]+1));
            return 0;
        }
    }

    g_tree_insert(ignore, GINT_TO_POINTER(idx[0]), GINT_TO_POINTER(idx[0]+1));
    g_tree_insert(safe, GINT_TO_POINTER(idx[1]), GINT_TO_POINTER(idx[1]+1));
    return 0;
}

gboolean same_point(gpointer key, gpointer value, gpointer data){
    args *comp_args;
    seg_3 *seg;
    comp_args = (args *)data;
    seg = comp_args->segs[GPOINTER_TO_INT(value)-1];
    if(seg->vert[0] == comp_args->comp || seg->vert[1] == comp_args->comp)
        *(comp_args->count) = *(comp_args->count) + 1;
    return false;
}

gboolean ig_tri_inner(gpointer key, gpointer value, gpointer data){

    d_args *args = (d_args *)data;
    seg_3 *outer, *inner;
    int idx = -1;

    if(!*(args->at_key)){
        if(args->key == key)
            *(args->at_key) = true;
        return false;
    }

    outer = args->segs[GPOINTER_TO_INT(args->key)];
    inner = args->segs[GPOINTER_TO_INT(key)];

    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
            if(outer->vert[i] == inner->vert[j])
                idx = segFromI(fmin(outer->idx[i+1 % 2], inner->idx[j+1 % 2]), fmax(outer->idx[i+1 % 2], inner->idx[j+1 % 2]), args->num_v);

    if(idx == -1 || GINT_TO_POINTER(idx) <= key || GINT_TO_POINTER(idx) <= args->key)
        return false;

    if(g_tree_lookup(args->tree, GINT_TO_POINTER(idx)))
        *(args->lc) = *(args->lc) + 1;


    return false;
}

gboolean ig_tri_outer(gpointer key, gpointer value, gpointer data){

    d_args *args = (d_args *)data;
    int local = 0;
    bool at_key = false;

    args->lc = &local;
    args->at_key = &at_key;
    args->key = key;

    g_tree_foreach(args->tree, ig_tri_inner, (gpointer)args);

    *(args->gc) = *(args->gc) + local;

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
    d_args *ig_tri_args;


    gsl_vector **vert;
    gsl_vector *x;
    seg_3 **seg;
    GTree *safe, *ignore;

    int currTri;
    tri_3 **tri;

    tri_share *seg_in_tri;
    int l;
    int num_t_use;
    // triangles per segment, ignored per segment
    int tps, ips;

    if(argc != 3){
        //n is number of vertices
        //file is csv of 3d vertices
        printf("USAGE: shapegen <n> <file>\n");
        exit(-1);
    }

    tri_ol = malloc(sizeof(args));
    ig_tri_args = malloc(sizeof(d_args));

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
            seg[iter]->vert[0] = vert[i];
            seg[iter]->vert[1] = vert[j];
            seg[iter]->idx[0] = i;
            seg[iter]->idx[1] = j;
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
                    printf("Seg %d (%d, %d) and seg %d (%d, %d)\n",
                            i, seg[i]->idx[0], seg[i]->idx[1],
                            j, seg[j]->idx[0], seg[j]->idx[1]);
                    resolve_overlap(safe, ignore, i, j);
                }
            }
            gsl_vector_free(x);
        }
    }

    printf("Safe: \n");
    g_tree_foreach(safe, print_g_tree, (gpointer)seg);
    printf("\n");
    printf("Ignore: \n");
    g_tree_foreach(ignore, print_g_tree, (gpointer)seg);
    printf("\n");


    //Each triangle needs 3 unique segments
    //However the entire search space isn't covered
    num_t = (num_v)*(num_v-1)*(num_v-2)/6;
    printf("num_t: %d\n", num_t);
    num_t -= g_tree_nnodes(ignore)*(num_v-2);
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

    count = 0;
    ig_tri_args->segs = seg;
    ig_tri_args->tree = ignore;
    ig_tri_args->gc = &count;
    ig_tri_args->num_v = num_v;
    g_tree_foreach(ignore, ig_tri_outer, (gpointer)ig_tri_args);
    num_t -= count;
    printf("Triangles in ignore: %d\n", count);


    printf("Number of triangles: %d\n", num_t);

    seg_in_tri = malloc(num_s*sizeof(tri_share));

    for(int i = 0; i < num_s; i++){
        // Each segment is shared by at most num_v-2 triangles
        seg_in_tri[i].tri = malloc((num_v-2)*sizeof(int));
        seg_in_tri[i].n = 0;
    }

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

                // TODO: Find a better way to do this

                printf("Gen tri %2d: (%d, %d, %d)\n", currTri, i, j, k);

                tri[currTri] = malloc(sizeof(tri_3));
                tri[currTri]->use = false;
                tri[currTri]->ignore = false;
                tri[currTri]->vert[0] = vert[i];
                tri[currTri]->idx[0] = i;
                tri[currTri]->vert[1] = vert[j];
                tri[currTri]->idx[1] = j;
                tri[currTri]->vert[2] = vert[k];
                tri[currTri]->idx[2] = k;

                l = segFromI(i, j, num_v);
                seg_in_tri[l].tri[seg_in_tri[l].n] = currTri;
                seg_in_tri[l].n++;

                l = segFromI(i, k, num_v);
                seg_in_tri[l].tri[seg_in_tri[l].n] = currTri;
                seg_in_tri[l].n++;

                l = segFromI(j, k, num_v);
                seg_in_tri[l].tri[seg_in_tri[l].n] = currTri;
                seg_in_tri[l].n++;

                currTri++;
            }
        }
    }

    if(currTri != num_t)
        exit(-1);

    for(int i = 0; i < num_t; i++)
        printf("Tri %2d: (%d, %d, %d)\n", i, tri[i]->idx[0], tri[i]->idx[1], tri[i]->idx[2]);

    for(int i = 0; i < num_s; i++){
        printf("Seg %2d (%d, %d): ", i, seg[i]->idx[0], seg[i]->idx[1]);
        for(int j = 0; j < seg_in_tri[i].n; j++){
            printf("%d (%d, %d, %d); ", seg_in_tri[i].tri[j],
                    tri[seg_in_tri[i].tri[j]]->idx[0],
                    tri[seg_in_tri[i].tri[j]]->idx[1],
                    tri[seg_in_tri[i].tri[j]]->idx[2]);
        }
        printf("\n");
    }

    for(int i = 0; i < num_s; i++){
        if(seg_in_tri[i].n != 2)
            continue;
        for(int j = 0; j < seg_in_tri[i].n; j++){
            tri[seg_in_tri[i].tri[j]]->use = true;
        }
    }

    num_t_use = 0;
    for(int i = 0; i < num_t; i++)
        if(tri[i]->use)
            num_t_use++;

    printf("%d/%d\n", num_t_use, (num_v-2)*2);
    while(num_t_use != (num_v-2)*2){
        // For each segment, make sure that only 2 triangles are using it
        for(int i = 0; i < num_s; i++){
            tps = 0;
            ips = 0;
            if(seg_in_tri[i].n == 2)
                continue;

            for(int j = 0; j < seg_in_tri[i].n; j++){
                if(tri[seg_in_tri[i].tri[j]]->use)
                    tps++;
                if(tri[seg_in_tri[i].tri[j]]->ignore)
                    ips++;
            }

//            if((tps != 2 && ips != (seg_in_tri[i].n - 2)) || (tps == 2 && ips == (seg_in_tri[i].n - 2)))
//                continue;

            if(tps < 2)
                for(int j = 0; j < seg_in_tri[i].n; j++)
                    if(!tri[seg_in_tri[i].tri[j]]->ignore)
                        tri[seg_in_tri[i].tri[j]]->use = true;

            if(ips < (seg_in_tri[i].n - 2))
                for(int j = 0; j < seg_in_tri[i].n; j++)
                    if(!tri[seg_in_tri[i].tri[j]]->use)
                        tri[seg_in_tri[i].tri[j]]->ignore = true;

            if(ips > seg_in_tri[i].n - 2){
                printf("fewer than 2 tri's/segment\n");
                for(int j = 0; j < seg_in_tri[i].n; j++){
                    //tri[seg_in_tri[i].tri[j]]->use = false;
                    tri[seg_in_tri[i].tri[j]]->ignore = false;
                }
            }
            if(tps > 2){
                printf("Greater than 2 tri's/segment\n");
                for(int j = 0; j < seg_in_tri[i].n; j++){
                    tri[seg_in_tri[i].tri[j]]->use = false;
                    //tri[seg_in_tri[i].tri[j]]->ignore = false;
                }
            }
        }

        num_t_use = 0;
        for(int i = 0; i < num_t; i++)
            if(tri[i]->use)
                num_t_use++;

    }

    printf("Num tri used: %d/%d\n", num_t_use, (num_v-2)*2);

    printf("Use these triangles: ");
    for(int i = 0; i < num_t; i++){
        if(tri[i]->use)
            printf("%d ", i);
    }
    printf("\n");

    //Memory management section

    for(int i = num_t-1; i >= 0; i--)
        free(tri[i]);
    free(tri);

    for(int i = 0; i < num_s; i++)
        free(seg_in_tri[i].tri);
    free(seg_in_tri);

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

    free(ig_tri_args);
    free(tri_ol);
    return(0);
}
