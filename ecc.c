#include "ecc.h"
#include <time.h>

ECC* ecc_new(Elliptic_Curve* curve, const Point* G) {
    ECC* ecc = malloc(sizeof(ECC));
    if (!ecc) return NULL;
    ecc->curve = curve;
    mpz_t zero;
    mpz_init(zero);
    mpz_set_str(zero, "0", 10);
    ecc->G = point_new(zero, zero);
    ecc->G.is_infinity = G->is_infinity;
    mpz_set(ecc->G.x, G->x);
    mpz_set(ecc->G.y, G->y);
    mpz_clear(zero);
    if (!elliptic_curve_is_on_curve(curve, &ecc->G)) {
        fprintf(stderr, "Point G is not on the curve\n");
        ecc_free(ecc);
        return NULL;
    }
    return ecc;
}

void ecc_free(ECC* ecc) {
    if (ecc) {
        point_free(&ecc->G);
        free(ecc);
    }
}

static void generate_random_private_key(mpz_t result, const mpz_t p) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    mpz_t one;
    mpz_init(one);
    mpz_set_str(one, "1", 10);
    mpz_sub(result, p, one);
    mpz_urandomm(result, state, result);
    mpz_add(result, result, one);

    mpz_clear(one);
    gmp_randclear(state);

}

Key ecc_generate_key(ECC* ecc) {
    Key key;
    mpz_init(key.private_key);
    do {
        generate_random_private_key(key.private_key, ecc->curve->p);
        key.public_key = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, key.private_key);
    } while (!elliptic_curve_is_on_curve(ecc->curve, &key.public_key));
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

    point_free(&temp);
    mpz_clear(random_k);
    return (Ciphertext){P1, P2};

}

Point ecc_decrypt(ECC* ecc, const Ciphertext* ciphertext, const mpz_t private_key) {
    mpz_t zero;
    mpz_init(zero);
    mpz_set_str(zero, "0", 10);

    Point P1 = point_new(zero, zero);
    P1.is_infinity = ciphertext->P1.is_infinity;
    mpz_set(P1.x, ciphertext->P1.x);
    mpz_set(P1.y, ciphertext->P1.y);

    Point P2 = point_new(zero, zero);
    P2.is_infinity = ciphertext->P2.is_infinity;
    mpz_set(P2.x, ciphertext->P2.x);
    mpz_set(P2.y, ciphertext->P2.y);

    mpz_clear(zero);

    mpz_t neg_y;
    mpz_init(neg_y);
    mpz_sub(neg_y, ecc->curve->p, P1.y);
    Point neg_P1 = point_new(P1.x, neg_y);
    neg_P1.is_infinity = P1.is_infinity;

    Point temp = elliptic_curve_scalar_multiply(ecc->curve, &neg_P1, private_key);
    Point result = elliptic_curve_add(ecc->curve, &P2, &temp);

    point_free(&neg_P1);
    point_free(&temp);
    point_free(&P1);
    point_free(&P2);
    mpz_clear(neg_y);

    return result;

}

void key_free(Key *key) {
    mpz_clear(key->private_key);
    point_free(&key->public_key);
}