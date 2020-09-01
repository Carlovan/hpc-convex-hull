/****************************************************************************
 *
 * mpi-convex-hull.c
 *
 * Compute the convex hull of a set of points in 2D, parallelized using MPI
 *
 * Copyright (C) 2020 Giulio Carlassare <giulio.carlassare(at)studio.unibo.it>
 * Last updated on 2020-08-28
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
 * wrapping"; è stato usato MPI per parallelizzare l'esecuzione.
 *   Le coordinate dei vertici dell'inviluppo con il minor numero di
 * punti sono stampate su standard output.  Si veda la specifica del
 * progetto sul sito del corso per una descrizione completa:
 *
 * http://moreno.marzolla.name/teaching/HPC/
 *
 * Per compilare:
 *
 * mpicc -D_XOPEN_SOURCE=600 -std=c99 -Wall -Wpedantic -O2 src/mpi-convex-hull.c -o build/mpi-convex-hull -lm
 *
 * (il flag -D_XOPEN_SOURCE=600 e' superfluo perche' viene settato
 * nell'header "hpc.h", ma definirlo tramite la riga di comando fa si'
 * che il programma compili correttamente anche se non si include
 * "hpc.h", o per errore non lo si include come primo file).
 *
 * Per eseguire il programma si puo' usare la riga di comando:
 *
 * mpirun -n ... build/mpi-convex-hull < ace.in > ace.hull
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

MPI_Datatype point_mpi_t;

/**
 * Returns true if b is better than a, considering cur as the
 * last point in the convex hull. If the points are collinear,
 * b is better than a if it is further from cur.
 */
bool better_point(const point_t cur, const point_t a, const point_t b) {
    int t = turn(cur, a, b);
    return t == LEFT || (t == COLLINEAR && fcmp(dist(cur, a) + dist(a, b), dist(cur, b)) == 0);
}

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 * pset only contains the partition use dby the current process.
 * totalPoints is the total number of points from the input.
 */
void convex_hull(const points_t *pset, points_t *hull, const int totalPoints, const int procCount, const int rank)
{
    hull->n = 0;
    /* There can be at most totalPoints points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(totalPoints * sizeof(*(hull->p))); assert(hull->p);
    
    const int n = pset->n;
    const point_t *p = pset->p;
    point_t *partialResults = (point_t*)malloc(sizeof(point_t) * procCount); assert(partialResults);
    
    // Identify the leftmost point of the local partition
    point_t leftmost = p[0];
    for (int i = 1; i<n; i++) {
        if (p[i].x < leftmost.x || (p[i].x == leftmost.x && p[i].y > leftmost.y)) {
            leftmost = p[i];
        }
    }
    MPI_Allgather(&leftmost, 1, point_mpi_t, partialResults, 1, point_mpi_t, MPI_COMM_WORLD);

    // Find the global leftmost point
    leftmost = partialResults[0];
    for (int i = 1; i<procCount; i++) {
        if (partialResults[i].x < leftmost.x || (partialResults[i].x == leftmost.x && partialResults[i].y > leftmost.y)) {
            leftmost = partialResults[i];
        }
    }

    point_t cur = leftmost;
    point_t next;
    
    do {
        // Add the current vertex to the hull
        assert(hull->n < totalPoints);
        hull->p[hull->n] = cur;
        hull->n++;
        
        // Search for the best vertex in local partition
        next = p[0];
        for (int j=1; j<n; j++) {
            if (better_point(cur, next, p[j])) {
                next = p[j];
            }
        }

        // Find global best, that is the next in che CH
        MPI_Allgather(&next, 1, point_mpi_t, partialResults, 1, point_mpi_t, MPI_COMM_WORLD);
        next = partialResults[0];
        for (int j=1; j<procCount; j++) {
            if (better_point(cur, next, partialResults[j])) {
                next = partialResults[j];
            }
        }
        assert(!points_eq(cur, next));
        cur = next;
    } while (!points_eq(leftmost, cur));
    
    /* Trim the excess space in the convex hull array */
    hull->p = (point_t*)realloc(hull->p, (hull->n) * sizeof(*(hull->p)));
    assert(hull->p); 
}

int main(int argc, char **argv )
{
    points_t pset, hull, allPoints;
    double tstart, elapsed;
    int procCount, rank;

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_contiguous(2, MPI_DOUBLE, &point_mpi_t);
    MPI_Type_commit(&point_mpi_t);

    // Read data and send to everyone
    if (rank == 0) {
        read_input(stdin, &allPoints);
    }
    tstart = hpc_gettime(); // Communication is part of the computation

    // Subdivide the points between all the processes
    MPI_Bcast(&allPoints.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int partSize = (allPoints.n + procCount - 1) / procCount;
    int* counts = (int*)malloc(sizeof(int) * procCount); assert(counts);
    int* starts = (int*)malloc(sizeof(int) * procCount); assert(starts);
    starts[0] = 0;
    counts[0] = min(allPoints.n, partSize);
    for (int i = 1; i < procCount; i++) {
        starts[i] = starts[i-1] + partSize;
        counts[i] = min(allPoints.n, starts[i] + partSize) - starts[i];
    }
    pset.p = (point_t*)malloc(sizeof(point_t) * counts[rank]); assert(pset.p);
    pset.n = counts[rank];

    // Use Scatterv because the last partition is probably smaller
    MPI_Scatterv(allPoints.p, counts, starts, point_mpi_t, pset.p, pset.n, point_mpi_t, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        free(allPoints.p);
    }

    // Find convex hull
    convex_hull(&pset, &hull, allPoints.n, procCount, rank);

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
