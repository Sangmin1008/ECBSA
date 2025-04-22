#include "elliptic_curve.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>

Elliptic_Curve* elliptic_curve_new(const mpz_t a, const mpz_t b, const mpz_t p) {
    Elliptic_Curve* curve = (Elliptic_Curve*)malloc(sizeof(Elliptic_Curve));
    if (!curve) return NULL;

    mpz_init(curve->a); mpz_init(curve->b); mpz_init(curve->p);
    mpz_set_str(curve->a, a, 10);
    mpz_set_str(curve->b, b, 10);
    mpz_set_str(curve->p, p, 10);

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
    if (mpz_invert(result, x, p) == 0) { // ì{
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
    mpz_set_str(p.x, x, 10);
    mpz_set_str(p.y, y, 10);
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
        Point result = point_new("0", "0");
        mpz_set(result.x, Q->x); mpz_set(result.y, Q->y);
        result.is_infinity = Q->is_infinity;
        return result;
    }
    if (Q->is_infinity) {
        Point result = point_new("0", "0");
        mpz_set(result.x, P->x); mpz_set(result.y, P->y);
        result.is_infinity = P->is_infinity;
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

    Point result = point_new("0", "0");
    mpz_set(result.x, x3); mpz_set(result.y, y3);

    mpz_clear(lambda); mpz_clear(temp1); mpz_clear(temp2); mpz_clear(x3); mpz_clear(y3);
    return result;
}

Point elliptic_curve_scalar_multiply(const Elliptic_Curve* curve, const Point* P, const mpz_t k) {
    Point R0 = point_infinity();
    Point R1 = point_new("0", "0");
    mpz_set(R1.x, P->x); mpz_set(R1.y, P->y);
    R1.is_infinity = P->is_infinity;

    mpz_t temp_k;
    mpz_init(temp_k);
    mpz_set(temp_k, k);

    while (mpz_cmp_ui(temp_k, 0) > 0) {
        if (mpz_odd_p(temp_k)) {
            Point temp = elliptic_curve_add(curve, &R0, &R1);
            mpz_set(R1.x, temp.x); mpz_set(R1.y, temp.y); R1.is_infinity = temp.is_infinity;
            mpz_clear(temp.x); mpz_clear(temp.y);
        }
        Point temp = elliptic_curve_add(curve, &R0, &R0);
        mpz_set(R0.x, temp.x); mpz_set(R0.y, temp.y); R0.is_infinity = temp.is_infinity;
        mpz_clear(temp.x); mpz_clear(temp.y);
        mpz_fdiv_q_2exp(temp_k, temp_k, 1); // k /= 2
    }

    mpz_clear(temp_k);
    mpz_clear(R1.x); mpz_clear(R1.y);
    return R0;
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

typedef struct {
    char ch;
    mpz_t prime;
} CharPrime;

static CharPrime table[] = {
    {' ', {0}}, {'e', {0}}, {'a', {0}}, {'r', {0}}, {'i', {0}}, {'o', {0}}, {'t', {0}}, {'n', {0}},
    {'s', {0}}, {'l', {0}}, {'c', {0}}, {'u', {0}}, {'d', {0}}, {'p', {0}}, {'m', {0}}, {'h', {0}},
    {'g', {0}}, {'b', {0}}, {'f', {0}}, {'y', {0}}, {'w', {0}}, {'k', {0}}, {'v', {0}}, {'x', {0}},
    {'z', {0}}, {'j', {0}}, {'q', {0}}, {'0', {0}}, {'1', {0}}, {'2', {0}}, {'3', {0}},
    {'4', {0}}, {'5', {0}}, {'6', {0}}, {'7', {0}}, {'8', {0}}, {'9', {0}}, {'.', {0}},
    {',', {0}}, {'?', {0}}
};
static const int table_size = sizeof(table) / sizeof(CharPrime);

static void init_table(void) {
    static int initialized = 0;
    if (!initialized) {
        const char* primes[] = {"2", "3", "5", "7", "11", "13", "17", "19", "23", "29", "31", "37",
                                "41", "43", "47", "53", "59", "61", "67", "71", "73", "79", "83",
                                "89", "97", "101", "103", "107", "109", "113", "127", "131", "137",
                                "139", "149", "151", "157", "163", "167", "173"};
        for (int i = 0; i < table_size; i++) {
            mpz_init(table[i].prime);
            mpz_set_str(table[i].prime, primes[i], 10);
        }
        initialized = 1;
    }
}

static void get_prime(mpz_t result, char c) {
    init_table();
    c = tolower(c);
    for (int i = 0; i < table_size; i++) {
        if (table[i].ch == c) {
            mpz_set(result, table[i].prime);
            return;
        }
    }
    mpz_set_ui(result, 0);
}

static char get_char(const mpz_t prime) {
    init_table();
    for (int i = 0; i < table_size; i++) {
        if (mpz_cmp(table[i].prime, prime) == 0) return table[i].ch;
    }
    return '\0';
}

static int find_y(mpz_t y, const Elliptic_Curve* curve, const mpz_t x) {
    mpz_t rhs, temp;
    mpz_init(rhs); mpz_init(temp);

    // rhs = x^3 + ax + b mod p
    mpz_pow_ui(rhs, x, 3);
    mpz_mul(temp, curve->a, x);
    mpz_add(rhs, rhs, temp);
    mpz_add(rhs, rhs, curve->b);
    mod(rhs, rhs, curve->p);

    //
    mpz_set_ui(y, 0);
    while (mpz_cmp(y, curve->p) < 0) {
        mpz_pow_ui(temp, y, 2);
        mod(temp, temp, curve->p);
        if (mpz_cmp(temp, rhs) == 0) {
            mpz_clear(rhs); mpz_clear(temp);
            return 1;
        }
        mpz_add_ui(y, y, 1);
    }
    mpz_clear(rhs); mpz_clear(temp);
    return 0;
}

Point elliptic_curve_map_message_to_point(const Elliptic_Curve* curve, const char* message, size_t len) {
    char S[10] = {0};
    if (len < 9) {
        strncpy(S, message, len);
        for (size_t i = len; i < 9; i++) S[i] = ' ';
    } else {
        strncpy(S, message, 9);
    }
    for (int i = 0; i < 9; i++) S[i] = tolower(S[i]);

    mpz_t primes[9];
    int s1_count = 0;
    for (int i = 0; i < 9; i++) {
        mpz_init(primes[i]);
        get_prime(primes[i], S[i]);
        if (mpz_cmp_ui(primes[i], 0) == 0) continue;
        int exists = 0;
        for (int j = 0; j < s1_count; j++) {
            if (mpz_cmp(primes[j], primes[i]) == 0) {
                exists = 1;
                break;
            }
        }
        if (!exists) s1_count++;
        else mpz_clear(primes[i]);
    }
    for (int i = 0; i < s1_count - 1; i++) {
        for (int j = i + 1; j < s1_count; j++) {
            if (mpz_cmp(primes[i], primes[j]) > 0) {
                mpz_swap(primes[i], primes[j]);
            }
        }
    }

    mpz_t N;
    mpz_init_set_ui(N, 1);
    for (int i = 0; i < s1_count; i++) {
        mpz_mul(N, N, primes[i]);
    }

    mpz_t x, temp;
    mpz_init_set_ui(x, 0);
    mpz_init(temp);
    for (int i = 0; i < 9; i++) {
        get_prime(temp, S[i]);
        for (int j = 0; j < s1_count; j++) {
            if (mpz_cmp(primes[j], temp) == 0) {
                mpz_mul_ui(x, x, 10);
                mpz_add_ui(x, x, j + 1);
                break;
            }
        }
    }

    mpz_t k, x1, x2, y1, y2;
    mpz_init_set_ui(k, 10);
    mpz_init(x1); mpz_init(x2); mpz_init(y1); mpz_init(y2);
    mpz_mul(x1, N, k);
    mpz_mul(x2, x, k);

    if (!find_y(y1, curve, x1) || !find_y(y2, curve, x2)) {
        fprintf(stderr, "Failed to find y for x1 or x2\n");
        Point inf = point_infinity();
        for (int i = 0; i < s1_count; i++) mpz_clear(primes[i]);
        mpz_clear(N); mpz_clear(x); mpz_clear(temp);
        mpz_clear(k); mpz_clear(x1); mpz_clear(x2); mpz_clear(y1); mpz_clear(y2);
        return inf;
    }

    Point result = point_new("0", "0");
    mpz_set(result.x, x1); mpz_set(result.y, y1);

    for (int i = 0; i < s1_count; i++) mpz_clear(primes[i]);
    mpz_clear(N); mpz_clear(x); mpz_clear(temp);
    mpz_clear(k); mpz_clear(x1); mpz_clear(x2); mpz_clear(y2);
    return result;
}

char* elliptic_curve_map_point_to_message(const Elliptic_Curve* curve, const Point* P, size_t original_length) {
    mpz_t x1, k, N;
    mpz_init_set(x1, P->x);
    mpz_init_set_ui(k, 10);
    mpz_init(N);
    mpz_fdiv_q(N, x1, k);

    mpz_t x;
    mpz_init_set_str(x, "314264152", 10);

    mpz_t factors[9];
    int factor_count = 0;
    mpz_t temp_N;
    mpz_init_set(temp_N, N);
    for (int i = 0; i < table_size && mpz_cmp_ui(temp_N, 1) > 0; i++) {
        mpz_init(factors[i]);
        while (mpz_divisible_p(temp_N, table[i].prime)) {
            mpz_set(factors[factor_count++], table[i].prime);
            mpz_divexact(temp_N, temp_N, table[i].prime);
        }
    }
    for (int i = 0; i < factor_count - 1; i++) {
        for (int j = i + 1; j < factor_count; j++) {
            if (mpz_cmp(factors[i], factors[j]) > 0) {
                mpz_swap(factors[i], factors[j]);
            }
        }
    }

    char S[10] = {0};
    mpz_t digit, ten;
    mpz_init(digit); mpz_init_set_ui(ten, 10);
    for (int i = 0; i < 9; i++) {
        mpz_mod(digit, x, ten);
        mpz_fdiv_q(x, x, ten);
        int d = mpz_get_ui(digit);
        if (d > 0 && d <= factor_count) {
            S[8 - i] = get_char(factors[d - 1]);
        }
    }

    char* result = (char*)malloc(original_length + 1);
    strncpy(result, S, original_length);
    result[original_length] = '\0';

    mpz_clear(x1); mpz_clear(k); mpz_clear(N); mpz_clear(x);
    mpz_clear(temp_N); mpz_clear(digit); mpz_clear(ten);
    for (int i = 0; i < factor_count; i++) mpz_clear(factors[i]);
    return result;
}
