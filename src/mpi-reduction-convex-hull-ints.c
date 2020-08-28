/****************************************************************************
 *
 * convex-hull.c
 *
 * Compute the convex hull of a set of points in 2D
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
#include <mpi.h>

#include "ioutils.h"
#include "data_types.h"
#include "mathutils.h"

void print_point(int rank, const char *text, const int p){
    printf("[%d] %s %d\n", rank, text, p);
}

// Check if b is better than a
bool better_point(const point_t cur, const point_t a, const point_t b) {
    int t = turn(cur, a, b);
    return t == LEFT || (t == COLLINEAR && fcmp(dist(cur, a) + dist(a, b), dist(cur, b)) == 0);
}

int gcur;
const point_t *gpoints;

void points_reduce_fn(void *invec_void, void *inoutvec_void, int *len, MPI_Datatype *dptr) {
    int *invec = (int*)invec_void;
    int *inoutvec = (int*)inoutvec_void;
    for (int i = 0; i < *len; i++) {
        if (better_point(gpoints[gcur], gpoints[inoutvec[i]], gpoints[invec[i]])) {
            inoutvec[i] = invec[i];
        }
    }
}
MPI_Op points_reduce_op;

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(const points_t *pset, points_t *hull, const int procCount, const int rank)
{
    const int n = pset->n;
    const point_t *p = pset->p;
    gpoints = p;

    /* There can be at most totalPoints points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->n = 0;
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);
    
    int *partialResults = (int*)malloc(sizeof(int) * procCount); assert(partialResults);

    int part_start = (n + procCount - 1) / procCount * rank;
    int part_end = min(n, (n + procCount - 1) / procCount * (rank + 1));
    
    // Identify the leftmost point of the local partition
    int leftmost = 0;
    for (int i = 1; i<n; i++) {
        if (p[i].x < p[leftmost].x || (p[i].x == p[leftmost].x && p[i].y > p[leftmost].y)) {
            leftmost = i;
        }
    }
    MPI_Allgather(&leftmost, 1, MPI_INT, partialResults, 1, MPI_INT, MPI_COMM_WORLD);

    // Find the global leftmost point
    leftmost = partialResults[0];
    for (int i = 1; i<procCount; i++) {
        if (p[partialResults[i]].x < p[leftmost].x || (p[partialResults[i]].x == p[leftmost].x && p[partialResults[i]].y > p[leftmost].y)) {
            leftmost = partialResults[i];
        }
    }
    print_point(rank, "leftmost", leftmost);
    int cur = leftmost;
    int next;
    
    do {
        /* Add the current vertex to the hull */
        assert(hull->n < n);
        hull->p[hull->n] = p[cur];
        hull->n++;
        
        // Search for the best vertex in local partition
        next = part_start;
        for (int j=part_start+1; j<part_end; j++) {
            if (better_point(p[cur], p[next], p[j])) {
                next = j;
            }
        }

        // Find global best, that is the next in che CH
        gcur = cur;
        int newNext;
        MPI_Allreduce(&next, &newNext, 1, MPI_INT, points_reduce_op, MPI_COMM_WORLD);
        next = newNext;
        assert(cur != next);
        cur = next;
    } while (leftmost != cur);
    
    /* Trim the excess space in the convex hull array */
    hull->p = (point_t*)realloc(hull->p, (hull->n) * sizeof(*(hull->p)));
    assert(hull->p); 
}

int main(int argc, char **argv )
{
    points_t pset, hull;
    double tstart, elapsed;
    int procCount, rank;

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Op_create(points_reduce_fn, true, &points_reduce_op);

    // Read data and send to everyone
    if (rank == 0) {
        read_input(stdin, &pset);
    }
    tstart = hpc_gettime(); // Communication is part of the computation

    MPI_Bcast(&pset.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        pset.p = (point_t*)malloc(sizeof(point_t) * pset.n); assert(pset.p);
    }
    MPI_Bcast(pset.p, 2*pset.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Find convex hull
    convex_hull(&pset, &hull, procCount, rank);

    elapsed = hpc_gettime() - tstart;

    if (rank == 0) {
        print_info(pset, hull, elapsed);
        write_hull(stdout, &hull);
    }
    free_pointset(&pset);
    free_pointset(&hull);
    MPI_Finalize();
    return EXIT_SUCCESS;    
}
