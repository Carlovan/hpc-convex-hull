#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "data_types.h"
#include <math.h>

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
 * Returns the dot product between the vectors p0--p1 and p1--p2
 */
double consecutive_dot_prod(const point_t p0, const point_t p1, const point_t p2) {
    const double x1 = p2.x - p1.x;
    const double y1 = p2.y - p1.y;
    const double x2 = p1.x - p0.x;
    const double y2 = p1.y - p0.y;
    return x1*x2 + y1*y2;
}

/**
 * Get the clockwise angle between the line p0p1 and the vector p1p2 
 *
 *         .
 *        . 
 *       .--+ (this angle) 
 *      .   |    
 *     .    V
 *    p1--------------p2
 *    /
 *   /
 *  /
 * p0
 *
 */
double cw_angle(const point_t p0, const point_t p1, const point_t p2)
{
    const double x1 = p2.x - p1.x;
    const double y1 = p2.y - p1.y;
    const double x2 = p1.x - p0.x;
    const double y2 = p1.y - p0.y;
    const double dot = x1*x2 + y1*y2;
    const double det = x1*y2 - y1*x2;
    const double result = atan2(det, dot);
    return (result >= 0 ? result : 2*M_PI + result);
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

#endif
