#include "elliptic_curve.h"


Elliptic_Curve* elliptic_curve_new(const mpz_t a, const mpz_t b, const mpz_t p) {
    Elliptic_Curve* curve = malloc(sizeof(Elliptic_Curve));
    if (!curve) return NULL;

    mpz_init(curve->a);
    mpz_set(curve->a, a);
    mpz_init(curve->b);
    mpz_set(curve->b, b);
    mpz_init(curve->p);
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
        mpz_clear(curve->a);
        mpz_clear(curve->b);
        mpz_clear(curve->p);
        free(curve);
    }
}

static void mod(mpz_t result, const mpz_t x, const mpz_t p) {
    mpz_mod(result, x, p);
}

static void mod_inv(mpz_t result, const mpz_t x, const mpz_t p) {
    if (!mpz_invert(result, x, p)) {
        fprintf(stderr, "No modular inverse\n");
        mpz_set_str(result, "0", 10);
    }
}

bool elliptic_curve_is_valid_curve(const Elliptic_Curve* curve) {
    mpz_t temp1, temp2, result;
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(result);

    mpz_pow_ui(temp1, curve->a, 3);
    mpz_mul_ui(temp1, temp1, 4);
    mpz_pow_ui(temp2, curve->b, 2);
    mpz_mul_ui(temp2, temp2, 27);
    mpz_add(result, temp1, temp2);
    mod(result, result, curve->p);

    bool valid = (mpz_cmp_ui(result, 0) != 0);

    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(result);
    return valid;

}

Point point_new(const mpz_t x, const mpz_t y) {
    Point p;
    mpz_init(p.x);
    mpz_set(p.x, x);
    mpz_init(p.y);
    mpz_set(p.y, y);
    p.is_infinity = false;
    return p;
}

Point point_infinity(void) {
    Point p;
    mpz_init(p.x);
    mpz_set_str(p.x, "0", 10);
    mpz_init(p.y);
    mpz_set_str(p.y, "0", 10);
    p.is_infinity = true;
    return p;
}

bool point_equals(const Point* p1, const Point* p2) {
    if (p1->is_infinity && p2->is_infinity) return true;
    if (!p1->is_infinity && !p2->is_infinity &&
        mpz_cmp(p1->x, p2->x) == 0 && mpz_cmp(p1->y, p2->y) == 0)
        return true;
    return false;
}

Point elliptic_curve_add(const Elliptic_Curve* curve, const Point* P, const Point* Q) {
    if (P->is_infinity) return point_new(Q->x, Q->y);
    if (Q->is_infinity) return point_new(P->x, P->y);
    if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) != 0) return point_infinity();

    mpz_t lambda, t1, t2, x3, y3, zero;
    mpz_init(lambda); mpz_init(t1); mpz_init(t2); mpz_init(x3); mpz_init(y3);
    mpz_init(zero);
    mpz_set_str(zero, "0", 10);

    if (point_equals(P, Q)) {
        mpz_pow_ui(t1, P->x, 2);
        mpz_mul_ui(t1, t1, 3);
        mpz_add(t1, t1, curve->a);
        mpz_mul_ui(t2, P->y, 2);
        mod_inv(t2, t2, curve->p);
        mpz_mul(lambda, t1, t2);
    } else {
        mpz_sub(t1, Q->y, P->y);
        mpz_sub(t2, Q->x, P->x);
        mod_inv(t2, t2, curve->p);
        mpz_mul(lambda, t1, t2);
    }
    mod(lambda, lambda, curve->p);

    mpz_pow_ui(x3, lambda, 2);
    mpz_sub(x3, x3, P->x);
    mpz_sub(x3, x3, Q->x);
    mod(x3, x3, curve->p);

    mpz_sub(y3, P->x, x3);
    mpz_mul(y3, lambda, y3);
    mpz_sub(y3, y3, P->y);
    mod(y3, y3, curve->p);

    Point R = point_new(x3, y3);

    mpz_clear(lambda); mpz_clear(t1); mpz_clear(t2);
    mpz_clear(x3); mpz_clear(y3); mpz_clear(zero);
    return R;

}

Point elliptic_curve_scalar_multiply(const Elliptic_Curve* curve, const Point* P, const mpz_t k) {
    Point R = point_infinity();
    Point Q = point_new(P->x, P->y);
    Q.is_infinity = P->is_infinity;

    mpz_t temp_k;
    mpz_init_set(temp_k, k);

    while (mpz_cmp_ui(temp_k, 0) > 0) {
        if (mpz_odd_p(temp_k)) {
            Point tmp = elliptic_curve_add(curve, &R, &Q);
            point_free(&R);
            R = tmp;
        }
        Point tmp2 = elliptic_curve_add(curve, &Q, &Q);
        point_free(&Q);
        Q = tmp2;
        mpz_fdiv_q_2exp(temp_k, temp_k, 1);
    }

    mpz_clear(temp_k);
    point_free(&Q);
    return R;

}

bool elliptic_curve_is_on_curve(const Elliptic_Curve* curve, const Point* P) {
    if (P->is_infinity) return true;

    mpz_t lhs, rhs, tmp;
    mpz_init(lhs); mpz_init(rhs); mpz_init(tmp);

    mpz_pow_ui(lhs, P->y, 2);
    mod(lhs, lhs, curve->p);

    mpz_pow_ui(rhs, P->x, 3);
    mpz_mul(tmp, curve->a, P->x);
    mpz_add(rhs, rhs, tmp);
    mpz_add(rhs, rhs, curve->b);
    mod(rhs, rhs, curve->p);

    bool on = (mpz_cmp(lhs, rhs) == 0);
    mpz_clear(lhs); mpz_clear(rhs); mpz_clear(tmp);
    return on;

}

void point_free(Point *p) {
    mpz_clear(p->x);
    mpz_clear(p->y);
}