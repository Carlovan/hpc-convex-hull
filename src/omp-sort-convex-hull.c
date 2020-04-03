/****************************************************************************
 *
 * omp-convex-hull.c
 *
 * Compute the convex hull of a set of points in 2D, parallelized with OpenMP
 *
 * Copyright (C) 2020 Giulio Carlassare <giulio.carlassare(at)studio.unibo.it>
 * Last updated on 2020-03-28
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************************
 *
 * Questo programma calcola l'inviluppo convesso (convex hull) di un
 * insieme di punti 2D letti da standard input usando l'algoritmo
 * "gift wrapping". Le coordinate dei vertici dell'inviluppo sono
 * stampate su standard output.  Per una descrizione completa del
 * problema si veda la specifica del progetto sul sito del corso:
 *
 * http://moreno.marzolla.name/teaching/HPC/
 *
 * Per compilare:
 *
 * gcc -D_XOPEN_SOURCE=600 -std=c99 -Wall -Wpedantic -O2 convex-hull.c -o convex-hull -lm
 *
 * (il flag -D_XOPEN_SOURCE=600 e' superfluo perche' viene settato
 * nell'header "hpc.h", ma definirlo tramite la riga di comando fa si'
 * che il programma compili correttamente anche se non si include
 * "hpc.h", o per errore non lo si include come primo file).
 *
 * Per eseguire il programma si puo' usare la riga di comando:
 *
 * ./convex-hull < ace.in > ace.hull
 * 
 * Per visualizzare graficamente i punti e l'inviluppo calcolato è
 * possibile usare lo script di gnuplot (http://www.gnuplot.info/)
 * incluso nella specifica del progetto:
 *
 * gnuplot -c plot-hull.gp ace.in ace.hull ace.png
 *
 ****************************************************************************/
#include "hpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "ioutils.h"
#include "data_types.h"
#include "mathutils.h"
#include <omp.h>
#include <string.h>

/**
 * Returns true if b is better than a
 */
inline bool better_point(const point_t cur, const point_t a, const point_t b) {
    int t = turn(cur, a, b);
    return t == LEFT || (t == COLLINEAR && consecutive_dot_prod(cur, a, b) > 0);
}

#pragma omp declare reduction ( leftmost_point : point_t* : \
        omp_out = omp_in->x < omp_out->x ? omp_in : omp_out ) initializer ( omp_priv = omp_orig )

void swap_points(point_t *a, point_t *b) {
    point_t tmp = *a;
    *a = *b;
    *b = tmp;
}

void swap_point_pointers(point_t **a, point_t ** b) {
    point_t *tmp = *a;
    *a = *b;
    *b = tmp;
}

point_t leftmostTmp;

/**
 * Returns < 0 if a < b
 */
int angle_cmp(const void* a, const void* b) {
    return turn(leftmostTmp, *(point_t*)b, *(point_t*)a);
}

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(points_t *pset, points_t *hull)
{
    const int n = pset->n;
    point_t *p = pset->p;

    hull->n = 0;
    /* There can be at most n points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);

    point_t *sortTmp = (point_t*)malloc(n * sizeof(*(pset->p))); assert(sortTmp);

    int num_threads;
    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();

    /* Identify the leftmost point p[leftmost] */
    point_t *leftmostP = &p[0];
    #pragma omp parallel for reduction(leftmost_point:leftmostP)
    for (int i = 1; i<n; i++) {
        if (p[i].x < leftmostP->x) {
            leftmostP = &p[i];
        }
    }
    leftmostTmp = *leftmostP; // Global, used fro comparison
    swap_points(leftmostP, &p[n-1]);

    // Sort points
    #pragma omp parallel
    {
        const int id = omp_get_thread_num();
        const int part_size = (n-1 + num_threads-1) / num_threads;
        const int start = id * part_size;
        int end = min((id + 1) * part_size, n-1);
        qsort(p+start, end-start, sizeof(*p), angle_cmp);

        #pragma omp barrier

        for (int len = 2; len/2 < num_threads; len *= 2) {
            if (id % len == 0) {
                end = min((id + len) * part_size, n-1);
                const int middle = min(end, (id + len/2) * part_size);
                int left = start, right = middle;
                int next = start;
                while (right < end || left < middle) {
                    assert(next < end);
                    if (right >= end || (left < middle && angle_cmp(&p[left], &p[right]) < 0)) {
                        sortTmp[next++] = p[left++];
                    } else {
                        sortTmp[next++] = p[right++];
                    }
                }
            }
            #pragma omp barrier
            #pragma omp single
            swap_point_pointers(&sortTmp, &p);
        }
    }
    p[n-1] = leftmostTmp;
    pset->p = p;
    free(sortTmp);

    // next[i] is the point after [i] in the hull
    int* next = (int*)malloc(n * sizeof(int)); assert(next);
    // prev[i] is the point before [i] in the hull
    int* prev = (int*)malloc(n * sizeof(int)); assert(prev);

    memset(next, -1, n * sizeof(int));

    int leftmost = n-1;

    #pragma omp parallel
    {
        const int id = omp_get_thread_num();
        const int part_size = (n + num_threads - 1) / num_threads;
        const int start = id * part_size;
        const int end = min((id + 1) * part_size, n-1);

        int cur = id == 0 ? leftmost : start;
        do {
            next[cur] = (cur+1) % n;
            for(int i = start; i <= end; i++) {
                if (better_point(p[cur], p[next[cur]], p[i])) {
                    next[cur] = i;
                }
            }
            cur = next[cur];
        } while(cur < end);
    }

    for(int i = 0; i < n; i++) {
        if (next[i] != -1)
            prev[next[i]] = i;
    }

    int cur = leftmost;
    do {
        bool removed = false;

        while(better_point(p[cur], p[next[cur]], p[next[next[cur]]])) {
            next[cur] = next[next[cur]];
            removed = true;
        }
        cur = removed ? prev[cur] : next[cur];
    } while(cur != leftmost);

    cur = leftmost;
    do {
        hull->p[hull->n] = p[cur];
        hull->n++;
        cur = next[cur];
    } while(cur != leftmost);

    free(next);
    free(prev);

    /* Trim the excess space in the convex hull array */
    hull->p = (point_t*)realloc(hull->p, (hull->n) * sizeof(*(hull->p)));
    assert(hull->p); 
}


int main( void )
{
    points_t pset, hull;
    double tstart, elapsed;

    read_input(stdin, &pset);
    tstart = hpc_gettime();
    convex_hull(&pset, &hull);
    elapsed = hpc_gettime() - tstart;
    print_info(pset, hull, elapsed, omp_get_max_threads());
    write_hull(stdout, &hull);
    free_pointset(&pset);
    free_pointset(&hull);
    return EXIT_SUCCESS;
}
