#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "data_types.h"
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int min(const int a, const int b) {
    return a < b ? a : b;
}

enum {
    LEFT = -1,
    COLLINEAR,
    RIGHT
};

/**
 * Return LEFT, RIGHT or COLLINEAR depending on the shape
 * of the vectors p0p1 and p1p2
 *
 * LEFT            RIGHT           COLLINEAR
 * 
 *  p2              p1----p2            p2
 *    \            /                   /
 *     \          /                   /
 *      p1       p0                  p1
 *     /                            /
 *    /                            /
 *  p0                            p0
 *
 * See Cormen, Leiserson, Rivest and Stein, "Introduction to Algorithms",
 * 3rd ed., MIT Press, 2009, Section 33.1 "Line-Segment properties"
 */
int turn(const point_t p0, const point_t p1, const point_t p2) {
    /*
      This function returns the correct result (COLLINEAR) also in the
      following cases:
      
      - p0==p1==p2
      - p0==p1
      - p1==p2
    */
    const double cross = (p1.x-p0.x)*(p2.y-p0.y) - (p2.x-p0.x)*(p1.y-p0.y);
    if (cross > 0.0) {
        return LEFT;
    } else {
        if (cross < 0.0) {
            return RIGHT;
        } else {
            return COLLINEAR;
        }
    }
}

/**
 * Compares two double values using a fixed similarity threshold.
 * Returns -1 if a < b, 0 if a == b and 1 if a > b.
 */
int fcmp(double a, double b) {
    const double t = 1e-10;
    if (fabs(a-b) < t) {
        return 0;
    }
    return (a > b) * 2 - 1;
}

double squared_dist(point_t p1, point_t p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return (dx*dx) + (dy*dy);
}

double dist(point_t p1, point_t p2) {
    return sqrt(squared_dist(p1, p2));
}

/**
 * Return true if the points are equal (compared usign fcmp).
 */
bool points_eq(const point_t a, const point_t b) {
    return fcmp(a.x, b.x) == 0 && fcmp(a.y, b.y) == 0;
}

#endif
