#include "elliptic_curve.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <C:\msys64\gmp-6.3.0\include\gmp.h>

Elliptic_Curve* elliptic_curve_new(const mpz_t a, const mpz_t b, const mpz_t p) {
    Elliptic_Curve* curve = (Elliptic_Curve*)malloc(sizeof(Elliptic_Curve));
    if (!curve) return NULL;

    mpz_init(curve->a); mpz_init(curve->b); mpz_init(curve->p);
    mpz_set(curve->a, a);
    mpz_set(curve->b, b);
    mpz_set(curve->p, p);

    if (!elliptic_curve_is_valid_curve(curve)) {
        fprintf(stderr, "Wrong Elliptic Curve\n");
        elliptic_curve_free(curve);
        return NULL;
    }
    return curve;
}

void elliptic_curve_free(Elliptic_Curve* curve) {
    if (curve) {
        mpz_clear(curve->a); mpz_clear(curve->b); mpz_clear(curve->p);
        free(curve);
    }
}

static void mod(mpz_t result, const mpz_t x, const mpz_t p) {
    mpz_mod(result, x, p);
}

static void mod_inv(mpz_t result, const mpz_t x, const mpz_t p) {
    if (mpz_invert(result, x, p) == 0) { // �{
        fprintf(stderr, "There is not mod_inv\n");
        mpz_set_ui(result, 0);
    }
}

bool elliptic_curve_is_valid_curve(const Elliptic_Curve* curve) {
    mpz_t temp1, temp2, result;
    mpz_init(temp1); mpz_init(temp2); mpz_init(result);

    // 4a^3
    mpz_pow_ui(temp1, curve->a, 3);
    mpz_mul_ui(temp1, temp1, 4);

    // 27b^2
    mpz_pow_ui(temp2, curve->b, 2);
    mpz_mul_ui(temp2, temp2, 27);

    // 4a^3 + 27b^2 mod p
    mpz_add(result, temp1, temp2);
    mod(result, result, curve->p);

    int is_valid = (mpz_cmp_ui(result, 0) != 0);
    mpz_clear(temp1); mpz_clear(temp2); mpz_clear(result);
    return is_valid;
}

Point point_new(const mpz_t x, const mpz_t y) {
    Point p;
    mpz_init(p.x); mpz_init(p.y);
    mpz_set(p.x, x);
    mpz_set(p.y, y);
    p.is_infinity = false;
    return p;
}

Point point_infinity(void) {
    Point p;
    mpz_init(p.x); mpz_init(p.y);
    mpz_set_ui(p.x, 0); mpz_set_ui(p.y, 0);
    p.is_infinity = true;
    return p;
}

bool point_equals(const Point* p1, const Point* p2) {
    return (p1->is_infinity && p2->is_infinity) ||
           (!p1->is_infinity && !p2->is_infinity &&
            mpz_cmp(p1->x, p2->x) == 0 && mpz_cmp(p1->y, p2->y) == 0);
}

Point elliptic_curve_add(const Elliptic_Curve* curve, const Point* P, const Point* Q) {
    if (P->is_infinity) {
        mpz_t zero; mpz_init_set_ui(zero, 0);
        Point result = point_new(zero, zero);
        mpz_set(result.x, Q->x); mpz_set(result.y, Q->y);
        result.is_infinity = Q->is_infinity;
        mpz_clear(zero);
        return result;
    }
    if (Q->is_infinity) {
        mpz_t zero; mpz_init_set_ui(zero, 0);
        Point result = point_new(zero, zero);
        mpz_set(result.x, P->x); mpz_set(result.y, P->y);
        result.is_infinity = P->is_infinity;
        mpz_clear(zero);
        return result;
    }

    if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) != 0) return point_infinity();

    mpz_t lambda, temp1, temp2, x3, y3;
    mpz_init(lambda); mpz_init(temp1); mpz_init(temp2); mpz_init(x3); mpz_init(y3);

    if (point_equals(P, Q)) {
        // lambda = (3x^2 + a) / (2y)
        mpz_pow_ui(temp1, P->x, 2); // x^2
        mpz_mul_ui(temp1, temp1, 3); // 3x^2
        mpz_add(temp1, temp1, curve->a); // 3x^2 + a
        mpz_mul_ui(temp2, P->y, 2); // 2y
        mod_inv(temp2, temp2, curve->p); // (2y)^-1
        mpz_mul(lambda, temp1, temp2);
        mod(lambda, lambda, curve->p);
    } else {
        // lambda = (y2 - y1) / (x2 - x1)
        mpz_sub(temp1, Q->y, P->y); // y2 - y1
        mpz_sub(temp2, Q->x, P->x); // x2 - x1
        mod_inv(temp2, temp2, curve->p); // (x2 - x1)^-1
        mpz_mul(lambda, temp1, temp2);
        mod(lambda, lambda, curve->p);
    }

    // x3 = lambda^2 - x1 - x2
    mpz_pow_ui(x3, lambda, 2);
    mpz_sub(x3, x3, P->x);
    mpz_sub(x3, x3, Q->x);
    mod(x3, x3, curve->p);

    // y3 = lambda * (x1 - x3) - y1
    mpz_sub(y3, P->x, x3);
    mpz_mul(y3, lambda, y3);
    mpz_sub(y3, y3, P->y);
    mod(y3, y3, curve->p);

    mpz_t zero; mpz_init_set_ui(zero, 0);
    Point result = point_new(zero, zero);
    mpz_set(result.x, x3); mpz_set(result.y, y3);

    mpz_clear(lambda); mpz_clear(temp1); mpz_clear(temp2); mpz_clear(x3); mpz_clear(y3);
    mpz_clear(zero);
    return result;
}

Point elliptic_curve_scalar_multiply(const Elliptic_Curve* curve, const Point* P, const mpz_t k) {
    Point R = point_infinity();
    Point Q = point_new(P->x, P->y);
    Q.is_infinity = P->is_infinity;

    mpz_t temp_k;
    mpz_init_set(temp_k, k);

    while (mpz_cmp_ui(temp_k, 0) > 0) {
        if (mpz_odd_p(temp_k)) {
            Point temp = elliptic_curve_add(curve, &R, &Q);
            mpz_clear(R.x); mpz_clear(R.y);
            R = temp;
        }
        Point temp = elliptic_curve_add(curve, &Q, &Q);
        mpz_clear(Q.x); mpz_clear(Q.y);
        Q = temp;
        mpz_fdiv_q_2exp(temp_k, temp_k, 1);
    }

    mpz_clear(temp_k);
    mpz_clear(Q.x); mpz_clear(Q.y);
    return R;
}

bool elliptic_curve_is_on_curve(const Elliptic_Curve* curve, const Point* P) {
    if (P->is_infinity) return true;

    mpz_t left, right, temp;
    mpz_init(left); mpz_init(right); mpz_init(temp);

    // left = y^2
    mpz_pow_ui(left, P->y, 2);
    mod(left, left, curve->p);

    // right = x^3 + ax + b
    mpz_pow_ui(right, P->x, 3);
    mpz_mul(temp, curve->a, P->x);
    mpz_add(right, right, temp);
    mpz_add(right, right, curve->b);
    mod(right, right, curve->p);

    int on_curve = (mpz_cmp(left, right) == 0);
    mpz_clear(left); mpz_clear(right); mpz_clear(temp);
    return on_curve;
}

void point_free(Point *p) {
    mpz_clear(p->x);
    mpz_clear(p->y);
}
