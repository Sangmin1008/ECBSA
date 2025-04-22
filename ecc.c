#include "ecc.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

ECC* ecc_new(Elliptic_Curve* curve, const Point* G) {
    ECC* ecc = (ECC*)malloc(sizeof(ECC));
    if (!ecc) return NULL;
    ecc->curve = curve;
    ecc->G = point_new("0", "0");
    mpz_set(ecc->G.x, G->x); mpz_set(ecc->G.y, G->y);
    ecc->G.is_infinity = G->is_infinity;
    if (!elliptic_curve_is_on_curve(curve, &ecc->G)) {
        fprintf(stderr, "Point G is not on the curve\n");
        ecc_free(ecc);
        return NULL;
    }
    return ecc;
}

void ecc_free(ECC* ecc) {
    if (ecc) {
        mpz_clear(ecc->G.x); mpz_clear(ecc->G.y);
        free(ecc);
    }
}

static void generate_random_private_key(mpz_t result, const mpz_t p) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, rand());

    mpz_t one;
    mpz_init_set_ui(one, 1);
    mpz_sub(result, p, one);
    mpz_urandomm(result, state, result);
    mpz_add(result, result, one);

    mpz_clear(one);
    gmp_randclear(state);
}

Key ecc_generate_key(ECC* ecc) {
    Key key;
    mpz_init(key.private_key);
    generate_random_private_key(key.private_key, ecc->curve->p);
    key.public_key = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, key.private_key);
    return key;
}

Ciphertext ecc_encrypt(ECC* ecc, const Point* message, const Point* public_key) {
    if (!elliptic_curve_is_on_curve(ecc->curve, public_key)) {
        fprintf(stderr, "Public Key is not on the curve\n");
        exit(1);
    }

    mpz_t random_k;
    mpz_init(random_k);
    Point P1;
    do {
        generate_random_private_key(random_k, ecc->curve->p);
        P1 = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, random_k);
    } while (mpz_cmp_ui(P1.x, 0) == 0 && mpz_cmp_ui(P1.y, 0) == 0);

    Point temp = elliptic_curve_scalar_multiply(ecc->curve, public_key, random_k);
    Point P2 = elliptic_curve_add(ecc->curve, message, &temp);

    mpz_clear(temp.x); mpz_clear(temp.y);
    mpz_clear(random_k);
    return (Ciphertext){P1, P2};
}

Point ecc_decrypt(ECC* ecc, const Ciphertext* ciphertext, const mpz_t private_key) {
    Point P1 = point_new("0", "0");
    mpz_set(P1.x, ciphertext->P1.x); mpz_set(P1.y, ciphertext->P1.y);
    P1.is_infinity = ciphertext->P1.is_infinity;

    Point P2 = point_new("0", "0");
    mpz_set(P2.x, ciphertext->P2.x); mpz_set(P2.y, ciphertext->P2.y);
    P2.is_infinity = ciphertext->P2.is_infinity;

    mpz_t neg_y;
    mpz_init(neg_y);
    mpz_sub(neg_y, ecc->curve->p, P1.y);
    Point neg_P1 = point_new("0", "0");
    mpz_set(neg_P1.x, P1.x); mpz_set(neg_P1.y, neg_y);
    neg_P1.is_infinity = P1.is_infinity;

    Point temp = elliptic_curve_scalar_multiply(ecc->curve, &neg_P1, private_key);
    Point result = elliptic_curve_add(ecc->curve, &P2, &temp);

    mpz_clear(neg_y);
    mpz_clear(P1.x); mpz_clear(P1.y);
    mpz_clear(P2.x); mpz_clear(P2.y);
    mpz_clear(neg_P1.x); mpz_clear(neg_P1.y);
    mpz_clear(temp.x); mpz_clear(temp.y);
    return result;
}

void key_free(Key *key) {
    mpz_clear(key->private_key);
    mpz_clear(key->public_key.x);
    mpz_clear(key->public_key.y);
}