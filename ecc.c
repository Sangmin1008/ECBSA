#include "ecc.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

ECC* ecc_new(Elliptic_Curve* curve, const Point* G) {
    ECC* ecc = (ECC*)malloc(sizeof(ECC));
    if (!ecc) return NULL;
    ecc->curve = curve;
    ecc->G = *G;
    if (!elliptic_curve_is_on_curve(curve, G)) {
        fprintf(stderr, "점 G가 곡선 위에 없음\n");
        free(ecc);
        return NULL;
    }
    return ecc;
}

void ecc_free(ECC* ecc) {
    if (ecc) free(ecc);
}

int64_t generate_random_private_key(int64_t p) {
    static int seeded = 0;
    if (!seeded) {
        srand(time(NULL));  // 시드 초기화 (한 번만 실행)
        seeded = 1;
    }
    return 1 + (rand() % (p - 1));
}

Key ecc_generate_key(ECC* ecc) {
    int64_t private_key = generate_random_private_key(ecc->curve->p);
    Point public_key = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, private_key);

    Key key = {private_key, public_key};
    return key;
}

Ciphertext ecc_encrypt(ECC* ecc, const Point* message, const Point* public_key) {
    if (!elliptic_curve_is_on_curve(ecc->curve, public_key)) {
        fprintf(stderr, "공개 키가 곡선 위에 없음\n");
        exit(1);
    }
    int64_t random_k;
    Point P1;
    do {
        random_k = generate_random_private_key(ecc->curve->p);
        P1 = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, random_k);
    } while (P1.x == 0 && P1.y == 0);
    Point temp = elliptic_curve_scalar_multiply(ecc->curve, public_key, random_k);
    Point P2 = elliptic_curve_add(ecc->curve, message, &temp);

    Ciphertext ciphertext = {P1, P2};
    return ciphertext;
}

Point ecc_decrypt(ECC* ecc, const Ciphertext* ciphertext, int64_t private_key) {
    Point P1 = ciphertext->P1;
    Point P2 = ciphertext->P2;
    Point neg_P1 = {P1.x, -P1.y, P1.is_infinity};
    Point temp = elliptic_curve_scalar_multiply(ecc->curve, &neg_P1, private_key);
    return elliptic_curve_add(ecc->curve, &P2, &temp);
}