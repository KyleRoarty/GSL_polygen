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

gboolean print_g_tree(gpointer key, gpointer value, gpointer data){
    if(value == NULL)
        return true;

    seg_3 **seg = (seg_3 **)data;
    int i = GPOINTER_TO_INT(value)-1;

    printf("%2d (%d, %d); ", i, seg[i]->idx[0], seg[i]->idx[1]);
    return false;
}

// x < y; sz = num points (cla)
// x, y are points that define a segment
int segFromI(int x, int y, int sz){
    return x*sz-x*(x+1)/2+y-(x+1);
}

// Used in a GLib function
// If b > a, search left half of tree
// if b < a, search right half of tree (I think...)
int sort_integers(gconstpointer a, gconstpointer b){
    int int_a, int_b;
    int_a = GPOINTER_TO_INT(a);
    int_b = GPOINTER_TO_INT(b);

    return int_a-int_b;
}

// Calculates slope of a segment, which is just the delta of the two
// segment ab, slope = b-a
int segment_slope_3(seg_3 *seg){
    int err;
    seg->slope = gsl_vector_calloc(V_DIM);

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

// Creates vertexes from csv file
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
// Checks if two segments intersect using QR decomposition from gsl
// The vectors are overconstrained, so QR decomp needs to be used
// it uses least-squares, so need to check if the residual significantly small
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

// Checks if the overlapping point between two segments
// is on the line segments
int overlap_in_bounds_3(seg_3 *seg1, seg_3 *seg2, gsl_vector *ol){

    seg_3 *segs[2];
    gsl_vector *sum;
    int ret;

    segs[0] = seg1;
    segs[1] = seg2;

    sum = gsl_vector_calloc(segs[0]->vert[0]->size);

    //  Check if ol is infinity, return -1 if it is
    for(int i = 0 ; i < ol->size; i++){
        if(gsl_isinf(gsl_vector_get(ol, i))){
            ret = GSL_FAILURE;
            goto cleanup;
        }
    }

    for(int j = 0; j < 2; j++){
        // slope*ol+start
        gsl_vector_memcpy(sum, segs[j]->slope);
        gsl_vector_scale(sum, gsl_vector_get(ol, 0));
        gsl_vector_add(sum, segs[j]->vert[0]);

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

// Determines which segment should be ignored in the mesh generation
// If both are in the 'safe' list, ignore one
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

    #ifdef DEBUG
    printf("s_a: %d, i_a: %d, s_b: %d, i_b: %d\n", safe_i[0], ignore_i[0], safe_i[1], ignore_i[1]);
    #endif

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
            g_tree_insert(safe, GINT_TO_POINTER(idx[(i+1) % 2]), GINT_TO_POINTER(idx[(i+1) % 2]+1));
            return 0;
        }
    }

    g_tree_insert(ignore, GINT_TO_POINTER(idx[0]), GINT_TO_POINTER(idx[0]+1));
    g_tree_insert(safe, GINT_TO_POINTER(idx[1]), GINT_TO_POINTER(idx[1]+1));
    return 0;
}

// Inner for loop; outer loop loops over vertexes, inner loop (this)
// loops over segments
// increments if the segment contains the vertex that the outer loop is on
gboolean same_vert(gpointer key, gpointer value, gpointer data){
    args *comp_args;
    seg_3 *seg;
    comp_args = (args *)data;
    seg = comp_args->segs[GPOINTER_TO_INT(value)-1];
    if(seg->vert[0] == comp_args->comp || seg->vert[1] == comp_args->comp)
        *(comp_args->count) = *(comp_args->count) + 1;
    return false;
}

// Inner loop for checking if the ignore list contains any triangles
// Checks all segments to see if it has a matching point with the outer segment
// If it does, it then checks if the third segment is in the ignore list
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
                idx = segFromI(fmin(outer->idx[(i+1) % 2], inner->idx[(j+1) % 2]), fmax(outer->idx[(i+1) % 2], inner->idx[(j+1) % 2]), args->num_v);

    if(idx == -1 || GINT_TO_POINTER(idx) <= key || GINT_TO_POINTER(idx) <= args->key)
        return false;

    if(g_tree_lookup(args->tree, GINT_TO_POINTER(idx)))
        *(args->lc) = *(args->lc) + 1;


    return false;
}

// Outer loop for checking if the ignore list contains any triangles
// Loops over segments, then calls the inner loop
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

    //Variables for g_tree_foreach arguments
    args *same_vert_args;
    d_args *ig_tri_args;


    // Array of vertices gotten from input file
    gsl_vector **vert;

    // Array of segments generated from vert
    seg_3 **seg;

    // Intersection point of two segments, if it exists
    gsl_vector *int_pt;

    // Trees used to hold segments that overlap. Ones in ignore are not used
    // in the generation of the shape
    GTree *safe, *ignore;

    // Triangles generated from vert
    tri_3 **tri;
    int currTri;

    // Holds list of triangles sharing a given segment (seg 0, 1, 2...)
    tri_share *tri_in_seg;

    // number of triangles used. Should be equal to 2num_v-4 upon finish
    int num_t_use;

    // triangles per segment, ignored triangles per segment
    int tps, ips;

    if(argc != 3){
        //n is number of vertices
        //file is csv of 3d vertices
        printf("USAGE: shapegen <n> <file.csv>\n");
        exit(-1);
    }
    //TODO: Add check for file

    same_vert_args = malloc(sizeof(args));
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
            for(int k = 0, vfs[2] = {i, j}; k < 2; k++){
                seg[iter]->vert[k] = vert[vfs[k]];
                seg[iter]->idx[k] = vfs[k];
            }
            if( (err = segment_slope_3(seg[iter])) ){
                printf("Error! From segment_slope_3\n");
            }
            iter++;
        }
    }

    for(int i = 0; i < num_s; i++){
        for(int j = i+1; j < num_s; j++){
            int_pt = gsl_vector_calloc(V_DIM-1);
            if( segment_intersect_3(seg[i], seg[j], int_pt) == GSL_SUCCESS ){
                if( overlap_in_bounds_3(seg[i], seg[j], int_pt) == GSL_SUCCESS ){
                    #ifdef DEBUG
                    printf("Seg %d (%d, %d) and seg %d (%d, %d)\n",
                            i, seg[i]->idx[0], seg[i]->idx[1],
                            j, seg[j]->idx[0], seg[j]->idx[1]);
                    #endif
                    resolve_overlap(safe, ignore, i, j);
                }
            }
            gsl_vector_free(int_pt);
        }
    }

    #ifdef DEBUG
    printf("Safe: \n");
    g_tree_foreach(safe, print_g_tree, (gpointer)seg);
    printf("\n");
    printf("Ignore: \n");
    g_tree_foreach(ignore, print_g_tree, (gpointer)seg);
    printf("\n");
    #endif


    //Each triangle needs 3 unique segments
    //However the entire search space isn't covered
    num_t = (num_v)*(num_v-1)*(num_v-2)/6;
    #ifdef DEBUG
    printf("num_t: %d\n", num_t);
    #endif
    num_t -= g_tree_nnodes(ignore)*(num_v-2);
    #ifdef DEBUG
    printf("num_t: %d\n", num_t);
    #endif

    //Set up args, call for loop
    same_vert_args->segs = seg;
    for(int i = 0; i < num_v; i++){
        count = 0;
        same_vert_args->count = &count;
        same_vert_args->comp = vert[i];
        g_tree_foreach(ignore, same_vert, (gpointer)same_vert_args);
        num_t += ((count)*(count-1)/2);
        #ifdef DEBUG
        printf("Count %d: %d\n", i, count);
        #endif
    }

    //Set up args, call for loop
    count = 0;
    ig_tri_args->segs = seg;
    ig_tri_args->tree = ignore;
    ig_tri_args->gc = &count;
    ig_tri_args->num_v = num_v;
    g_tree_foreach(ignore, ig_tri_outer, (gpointer)ig_tri_args);
    num_t -= count;

    #ifdef DEBUG
    printf("Triangles in ignore: %d\n", count);
    printf("Number of triangles: %d\n", num_t);
    #endif

    tri_in_seg = malloc(num_s*sizeof(tri_share));

    for(int i = 0; i < num_s; i++){
        // Each segment is shared by at most num_v-2 triangles
        tri_in_seg[i].tri = malloc((num_v-2)*sizeof(int));
        tri_in_seg[i].n = 0;
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

                #ifdef DEBUG
                printf("Gen tri %2d: (%d, %d, %d)\n", currTri, i, j, k);
                #endif

                tri[currTri] = malloc(sizeof(tri_3));
                tri[currTri]->use = false;
                tri[currTri]->ignore = false;
                for(int l = 0, m = 0, vft[3] = {i, j, k}; l < 3; l++){
                    tri[currTri]->vert[l] = vert[vft[l]];
                    tri[currTri]->idx[l] = vft[l];
                    m = segFromI(fmin(vft[l], vft[(l+1) % 3]), fmax(vft[l], vft[(l+1) % 3]), num_v);
                    tri_in_seg[m].tri[tri_in_seg[m].n] = currTri;
                    tri_in_seg[m].n++;
                }

                currTri++;
            }
        }
    }

    if(currTri != num_t)
        exit(-1);


    #ifdef DEBUG
    for(int i = 0; i < num_t; i++)
        printf("Tri %2d: (%d, %d, %d)\n", i, tri[i]->idx[0], tri[i]->idx[1], tri[i]->idx[2]);

    for(int i = 0; i < num_s; i++){
        printf("Seg %2d (%d, %d): ", i, seg[i]->idx[0], seg[i]->idx[1]);
        for(int j = 0; j < tri_in_seg[i].n; j++){
            printf("%d (%d, %d, %d); ", tri_in_seg[i].tri[j],
                    tri[tri_in_seg[i].tri[j]]->idx[0],
                    tri[tri_in_seg[i].tri[j]]->idx[1],
                    tri[tri_in_seg[i].tri[j]]->idx[2]);
        }
        printf("\n");
    }
    #endif

    for(int i = 0; i < num_s; i++){
        if(tri_in_seg[i].n != 2)
            continue;
        for(int j = 0; j < tri_in_seg[i].n; j++){
            tri[tri_in_seg[i].tri[j]]->use = true;
        }
    }

    num_t_use = 0;
    for(int i = 0; i < num_t; i++)
        if(tri[i]->use)
            num_t_use++;

    #ifdef DEBUG
    printf("%d/%d\n", num_t_use, (num_v-2)*2);
    #endif

    while(num_t_use != (num_v-2)*2){
        // For each segment, make sure that only 2 triangles are using it
        for(int i = 0; i < num_s; i++){
            tps = 0;
            ips = 0;
            if(tri_in_seg[i].n == 2)
                continue;

            for(int j = 0; j < tri_in_seg[i].n; j++){
                if(tri[tri_in_seg[i].tri[j]]->use)
                    tps++;
                if(tri[tri_in_seg[i].tri[j]]->ignore)
                    ips++;
            }

            if(tps < 2)
                for(int j = 0; j < tri_in_seg[i].n; j++)
                    if(!tri[tri_in_seg[i].tri[j]]->ignore)
                        tri[tri_in_seg[i].tri[j]]->use = true;

            if(ips < (tri_in_seg[i].n - 2))
                for(int j = 0; j < tri_in_seg[i].n; j++)
                    if(!tri[tri_in_seg[i].tri[j]]->use)
                        tri[tri_in_seg[i].tri[j]]->ignore = true;


            // These two conditions are bad, so just unflag everything and
            // run through again
            if(ips > tri_in_seg[i].n - 2){
                #ifdef DEBUG
                printf("fewer than 2 tri's/segment\n");
                #endif
                for(int j = 0; j < tri_in_seg[i].n; j++){
                    tri[tri_in_seg[i].tri[j]]->ignore = false;
                }
            }
            if(tps > 2){
                #ifdef DEBUG
                printf("Greater than 2 tri's/segment\n");
                #endif
                for(int j = 0; j < tri_in_seg[i].n; j++){
                    tri[tri_in_seg[i].tri[j]]->use = false;
                }
            }
        }

        num_t_use = 0;
        for(int i = 0; i < num_t; i++)
            if(tri[i]->use)
                num_t_use++;

    }

    printf("Use these triangles: ");
    for(int i = 0; i < num_t; i++){
        if(tri[i]->use)
            printf("%2d ", i);
    }
    printf("\n");

    //Memory management section

    for(int i = num_t-1; i >= 0; i--)
        free(tri[i]);
    free(tri);

    for(int i = num_s-1; i >= 0; i--)
        free(tri_in_seg[i].tri);
    free(tri_in_seg);

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
    free(same_vert_args);
    return(0);
}
