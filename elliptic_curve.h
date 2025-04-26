#ifndef ECC_CAPSTONE_ELLIPTIC_CURVE_H
#define ECC_CAPSTONE_ELLIPTIC_CURVE_H

#include <stdbool.h>
#include <stddef.h>
#include <C:\msys64\gmp-6.3.0\include\gmp.h>

typedef struct {
    mpz_t a, b, p;
} Elliptic_Curve;

typedef struct {
    mpz_t x, y;
    bool is_infinity;
} Point;

Elliptic_Curve* elliptic_curve_new(const mpz_t a, const mpz_t b, const mpz_t p);
void elliptic_curve_free(Elliptic_Curve* curve);
bool elliptic_curve_is_valid_curve(const Elliptic_Curve* curve);

Point point_new(const mpz_t x, const mpz_t y);
Point point_infinity(void);
bool point_equals(const Point* p1, const Point* p2);

Point elliptic_curve_add(const Elliptic_Curve* curve, const Point* P, const Point* Q);
Point elliptic_curve_scalar_multiply(const Elliptic_Curve* curve, const Point* P, const mpz_t k);
bool elliptic_curve_is_on_curve(const Elliptic_Curve* curve, const Point* P);

void point_free(Point *p);
#endif // ECC_CAPSTONE_ELLIPTIC_CURVE_H
