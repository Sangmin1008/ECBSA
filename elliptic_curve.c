#include "elliptic_curve.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

Elliptic_Curve* elliptic_curve_new(int64_t a, int64_t b, int64_t p) {
    Elliptic_Curve* curve = (Elliptic_Curve*)malloc(sizeof(Elliptic_Curve));
    if (!curve) return NULL;
    curve->a = a;
    curve->b = b;
    curve->p = p;
    if (!elliptic_curve_is_valid_curve(curve)) {
        fprintf(stderr, "잘못된 타원 곡선\n");
        free(curve);
        return NULL;
    }
    return curve;
}

void elliptic_curve_free(Elliptic_Curve* curve) {
    if (curve) free(curve);
}

int64_t mod(int64_t x, int64_t p) {
    return ((x % p) + p) % p;
}

int64_t mod_inv(int64_t x, int64_t p) {
    int64_t res = 1, y = p - 2;
    while (y) {
        if (y & 1) res = (res * x) % p;
        x = (x * x) % p;
        y >>= 1;
    }
    return res;
}

bool elliptic_curve_is_valid_curve(const Elliptic_Curve* curve) {
    return mod(4 * curve->a * curve->a * curve->a + 27 * curve->b * curve->b, curve->p) != 0;
}

Point point_new(int64_t x, int64_t y) {
    Point p = { x, y, false };
    return p;
}

Point point_infinity(void) {
    Point p = { 0, 0, true };
    return p;
}

bool point_equals(const Point* p1, const Point* p2) {
    return (p1->is_infinity && p2->is_infinity) ||
           (!p1->is_infinity && !p2->is_infinity && p1->x == p2->x && p1->y == p2->y);
}

Point elliptic_curve_add(const Elliptic_Curve* curve, const Point* P, const Point* Q) {
    if (P->is_infinity) return *Q;
    if (Q->is_infinity) return *P;

    if (P->x == Q->x && P->y != Q->y) return point_infinity();

    int64_t lambda;
    if (point_equals(P, Q))
        lambda = mod((3 * P->x * P->x + curve->a) * mod_inv(2 * P->y, curve->p), curve->p);
    else
        lambda = mod((Q->y - P->y) * mod_inv(Q->x - P->x, curve->p), curve->p);

    int64_t x3 = mod(mod(lambda * lambda, curve->p) - P->x - Q->x, curve->p);
    int64_t y3 = mod(mod(lambda * (P->x - x3), curve->p) - P->y, curve->p);

    return point_new(x3, y3);
}

Point elliptic_curve_scalar_multiply(const Elliptic_Curve* curve, const Point* P, int64_t k) {
    Point R0 = point_infinity();
    Point R1 = *P;

    for (int i = 63; i >= 0; --i) {
        if (k & (1LL << i)) {
            R1 = elliptic_curve_add(curve, &R0, &R1);
            R0 = elliptic_curve_add(curve, &R0, &R0);
        }
        else {
            R0 = elliptic_curve_add(curve, &R0, &R1);
            R1 = elliptic_curve_add(curve, &R1, &R1);
        }
    }
    return R0;
}

bool elliptic_curve_is_on_curve(const Elliptic_Curve* curve, const Point* P) {
    if (P->is_infinity) return true;
    return mod(P->y * P->y, curve->p) == mod(P->x * P->x * P->x + curve->a * P->x + curve->b, curve->p);
}