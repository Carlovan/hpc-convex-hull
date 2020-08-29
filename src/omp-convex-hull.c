/****************************************************************************
 *
 * convex-hull.c
 *
 * Compute the convex hull of a set of points in 2D
 *
 * Copyright (C) 2020 Giulio Carlassare <giulio.carlassare(at)studio.unibo.it>
 * Last updated on 2020-05-05
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
 * insieme di punti 2D letti da standard input con l'algoritmo "gift
 * wrapping";  è stato usato OpenMP per parallelizzare la soluzione.
 *   Le coordinate dei vertici dell'inviluppo con il minor numero di
 * punti sono stampate su standard output.  Si veda la specifica del
 * progetto sul sito del corso per una descrizione completa:
 *
 * http://moreno.marzolla.name/teaching/HPC/
 *
 * Per compilare:
 *
 * gcc -D_XOPEN_SOURCE=600 -std=c99 -Wall -Wpedantic -O2 omp-convex-hull.c -o omp-convex-hull -lm
 *
 * (il flag -D_XOPEN_SOURCE=600 e' superfluo perche' viene settato
 * nell'header "hpc.h", ma definirlo tramite la riga di comando fa si'
 * che il programma compili correttamente anche se non si include
 * "hpc.h", o per errore non lo si include come primo file).
 *
 * Per eseguire il programma si puo' usare la riga di comando:
 *
 * ./omp-convex-hull < ace.in > ace.hull
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
#include <string.h>
#include <stdbool.h>

#include "ioutils.h"
#include "data_types.h"
#include "mathutils.h"

#pragma omp declare reduction (leftmost_point : point_t : omp_out = fcmp(omp_out.x, omp_in.x) < 0 ? omp_out : omp_in) \
    initializer ( omp_priv = omp_orig )

/* If every thread works on less than `copyThreshold` points,
 * then they are not copied */
const int copyThreshold = 2e9;

bool better_point(const point_t cur, const point_t a, const point_t b) {
    int t = turn(cur, a, b);
    return t == LEFT || (t == COLLINEAR && consecutive_dot_prod(cur, a, b) > 0);
}

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(const points_t *pset, points_t *hull) {
    const int n = pset->n;

    hull->n = 0;
    /* There can be at most n points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);
    
    /* Identify the leftmost point p[leftmost] */
    point_t leftmost = pset->p[0];
    #pragma omp parallel for default(none) shared(pset) firstprivate(n) reduction(leftmost_point:leftmost)
    for (int i = 1; i<n; i++) {
        if (pset->p[i].x < leftmost.x) {
            leftmost = pset->p[i];
        }
    }
    point_t cur = leftmost;

    /* Retrieve number of threads to allocate enough memory
       for reduced values of individual threads */
    int num_threads;
    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();

    point_t* results = (point_t*)malloc(num_threads * sizeof(point_t));
    assert(results);
    
    #pragma omp parallel default(none) firstprivate(n, num_threads) shared(pset, cur, hull, leftmost, results)
    {
        /* Points are equally subdivided between threads */
        const int id = omp_get_thread_num();
        const int part_size = (n + num_threads - 1) / num_threads;
        const int start = part_size * id;
        const int end = min(n, start + part_size);
        const int memsize = (end-start)*sizeof(point_t);

        /* If there are enough points, every thread copies the data it will work with
           in its own local memory (assuming there are multiple NUMA nodes). */
        point_t *p;
        bool copied = false;
        if (part_size > copyThreshold) {
            p = (point_t*)malloc(memsize); assert(p);
            memcpy(p, pset->p + start, memsize);
            copied = true;
        } else {
            p = &pset->p[start];
        }

        /* Actual Gift Wrapping algorithm */
        do {
            /* Add the last found point to the convex hull */
            #pragma omp single
            {
                assert(hull->n < n);
                hull->p[hull->n] = cur;
                hull->n++;
            }

            /* Search for the next vertex; every thread works on its partition */
            point_t next = p[0];
            for (int j=start+1; j<end; j++) {
                if(better_point(cur, next, p[j-start])) {
                    next = p[j-start];
                }
            }
            /* Save locally reduced value in common array */
            results[id] = next;

            /* Serially reduce the results of individual threads */
            #pragma omp barrier
            #pragma omp single
            {
                next = results[0];
                for(int j = 1; j < num_threads; j++) {
                    if (better_point(cur, next, results[j])) {
                        next = results[j];
                    }
                }
                assert(!points_eq(cur, next));
                cur = next;
            }
        } while (!points_eq(cur, leftmost));

        /* If there were few points we didn't copy */
        if (copied) {
            free(p);
        }
    }
    free(results);
    
    /* Trim the excess space in the convex hull array */
    hull->p = (point_t*)realloc(hull->p, (hull->n) * sizeof(*(hull->p)));
    assert(hull->p); 
}

int main( void ) {
    points_t pset, hull;
    double tstart, elapsed;
    
    read_input(stdin, &pset);
    tstart = hpc_gettime();
    convex_hull(&pset, &hull);
    elapsed = hpc_gettime() - tstart;
    print_info(pset, hull, elapsed);
    write_hull(stdout, &hull);
    free_pointset(&pset);
    free_pointset(&hull);
    return EXIT_SUCCESS;    
}
