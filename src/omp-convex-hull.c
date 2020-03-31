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

/**
 * Returns true if b is better than a
 */
inline bool better_point(const point_t cur, const point_t a, const point_t b) {
    int t = turn(cur, a, b);
    return t == LEFT || (t == COLLINEAR && consecutive_dot_prod(cur, a, b) > 0);
}

typedef struct {
    const point_t *point;
    const point_t * const cur;
} reduction_value_t;

#pragma omp declare reduction ( best_point : reduction_value_t : \
        omp_out.point = better_point(*omp_out.cur, *omp_out.point, *omp_in.point) \
            ? omp_in.point : omp_out.point )\
            initializer ( omp_priv = omp_orig )

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(const points_t *pset, points_t *hull)
{
    const int n = pset->n;
    const point_t *p = pset->p;

    hull->n = 0;
    /* There can be at most n points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);
    
    /* Identify the leftmost point p[leftmost] */
    int leftmost = 0;
    for (int i = 1; i<n; i++) {
        if (p[i].x < p[leftmost].x) {
            leftmost = i;
        }
    }
    int cur = leftmost;
    
    /* Main loop of the Gift Wrapping algorithm. This is where most of
       the time is spent; therefore, this is the block of code that
       must be parallelized. */
    do {
        /* Add the current vertex to the hull */
        assert(hull->n < n);
        hull->p[hull->n] = p[cur];
        hull->n++;
        
        /* Search for the next vertex */
        reduction_value_t next = {&p[(cur + 1) % n], &p[cur]};
        #pragma omp parallel for reduction(best_point:next)
        for (int j=0; j<n; j++) {
            if (better_point(p[cur], *next.point, p[j])) {
                next.point = &p[j];
            }
        }
        int nextI = next.point - p;
        assert(cur != nextI);
        cur = nextI;
    } while (cur != leftmost);
    
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
