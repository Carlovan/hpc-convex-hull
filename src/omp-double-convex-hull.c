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

/**
 * Checks if the point p is above the line between left and right
 */
bool is_above_line(const point_t left, const point_t right, const point_t p) {
    return (right.y - left.y) / (right.x - left.x) * (p.x - left.x) < (p.y - left.y);
}

typedef struct {
    const point_t *point;
    const point_t * const cur;
} reduction_value_t;

void swap_points(point_t *a, point_t *b) {
    point_t tmp = *a;
    *a = *b;
    *b = tmp;
}

#pragma omp declare reduction ( best_point : reduction_value_t : \
        omp_out.point = better_point(*omp_out.cur, *omp_out.point, *omp_in.point) \
            ? omp_in.point : omp_out.point )\
            initializer ( omp_priv = omp_orig )

#pragma omp declare reduction ( leftmost_point : point_t* : \
        omp_out = omp_in->x < omp_out->x ? omp_in : omp_out ) initializer ( omp_priv = omp_orig )
#pragma omp declare reduction ( rightmost_point : point_t* : \
        omp_out = omp_in->x > omp_out->x ? omp_in : omp_out ) initializer ( omp_priv = omp_orig )

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(const points_t *pset, points_t *hull)
{
    const int n = pset->n;
    point_t *p = pset->p;

    hull->n = 0;
    /* There can be at most n points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);

    points_t rhull;
    rhull.n = 0;
    rhull.p = (point_t*)malloc(n * sizeof(*(rhull.p))); assert(rhull.p);
    
    /* Identify the leftmost point p[leftmost] */
    point_t *leftmostP = &p[0];
    point_t *rightmostP = &p[0];
    #pragma omp parallel for reduction(leftmost_point:leftmostP) reduction(rightmost_point:rightmostP)
    for (int i = 1; i<n; i++) {
        if (p[i].x < leftmostP->x) {
            leftmostP = &p[i];
        }
        if (p[i].x > rightmostP->x) {
            rightmostP = &p[i];
        }
    }
    int leftmost = n-1;
    int rightmost = 0;
    swap_points(leftmostP, &p[leftmost]);
    swap_points(rightmostP, &p[rightmost]);

    int leftIndex = 1, rightIndex = n-2;
    while(leftIndex < rightIndex) {
        bool isAboveL = is_above_line(p[leftmost], p[rightmost], p[leftIndex]);
        bool isAboveR = is_above_line(p[leftmost], p[rightmost], p[rightIndex]);
        
        if (!isAboveL && isAboveR) {
            swap_points(&p[leftIndex], &p[rightIndex]);
        }
        if (isAboveL) {
            leftIndex++;
        } else if (!isAboveR) { // Ensure only one is updated in an iteration
            rightIndex--;
        }
    }

    for(int from_left = 0; from_left < 2; from_left++) {
        int cur, last;
        int start, end;
        points_t *my_hull;
        if (from_left) {
            cur = leftmost;
            last = rightmost;
            start = 0;
            end = leftIndex;
            my_hull = hull;
        } else {
            cur = rightmost;
            last = leftmost;
            start = leftIndex;
            end = n;
            my_hull = &rhull;
        }

        /* Main loop of the Gift Wrapping algorithm. This is where most of
           the time is spent; therefore, this is the block of code that
           must be parallelized. */
        do {
            /* Add the current vertex to the hull */
            my_hull->p[my_hull->n] = p[cur];
            my_hull->n++;
            
            /* Search for the next vertex */
            reduction_value_t next = {&p[last], &p[cur]};
            #pragma omp parallel for reduction(best_point:next)
            for (int j=start; j<end; j++) {
                if (better_point(p[cur], *next.point, p[j])) {
                    next.point = &p[j];
                }
            }
            int nextI = next.point - p;
            assert(cur != nextI);
            cur = nextI;
        } while (cur != last);
    }

    // Add all points from the right hull to the left one
    for (int i = 0; i < rhull.n; i++) {
        hull->p[hull->n] = rhull.p[i];
        hull->n++;
    }
    
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
