#ifndef _DATA_TYPES_H_
#define _DATA_TYPES_H_

#include <stdlib.h>
#include <stdbool.h>

/* A single point */
typedef struct {
    double x, y;
} point_t;

/* An array of n points */
typedef struct {
    int n;      /* number of points     */
    point_t *p; /* array of points      */
} points_t;

bool point_equals(const point_t a, const point_t b) {
	return a.x == b.x && a.y == b.y;
}

/**
 * Free the memory allocated by structure pset.
 */
void free_pointset( points_t *pset )
{
    pset->n = 0;
    free(pset->p);
    pset->p = NULL;
}

#endif
