#ifndef ECC_CAPSTONE_ECC_H
#define ECC_CAPSTONE_ECC_H

#include "elliptic_curve.h"

typedef struct {
    Elliptic_Curve* curve;
    Point G;
} ECC;

typedef struct {
    mpz_t private_key;
    Point public_key;
} Key;

typedef struct {
    Point P1;
    Point P2;
} Ciphertext;

ECC* ecc_new(Elliptic_Curve* curve, const Point* G);
void ecc_free(ECC* ecc);
void generate_random_private_key(mpz_t result, const mpz_t p);
Key ecc_generate_key(ECC* ecc);
Ciphertext ecc_encrypt(ECC* ecc, const Point* message, const Point* public_key);
Point ecc_decrypt(ECC* ecc, const Ciphertext* ciphertext, const mpz_t private_key); // private_key도 mpz_t로
void key_free(Key *key);
#endif // ECC_CAPSTONE_ECC_H
