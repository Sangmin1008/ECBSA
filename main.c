#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <stdlib.h>
#include "elliptic_curve.h"
#include "ecc.h"

int main() {
    srand(time(NULL));
    mpz_t prime, a, b, x, y, mx, my;

    mpz_init(prime);
    mpz_set_str(prime, "115792089237316195423570985008687907853269984665640564039457584007908834671663", 10);
    mpz_init(a);
    mpz_set_str(a, "0", 10);
    mpz_init(b);
    mpz_set_str(b, "7", 10);

    mpz_init(x);
    mpz_set_str(x, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_init(y);
    mpz_set_str(y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Elliptic_Curve* curve = elliptic_curve_new(a, b, prime);
    Point G = point_new(x, y);
    mpz_mod(G.x, G.x, curve->p);
    mpz_mod(G.y, G.y, curve->p);

    ECC* ecc = ecc_new(curve, &G);
    if (!ecc) {
        fprintf(stderr, "Failed to create ECC\n");
        exit(1);
    }

    Key key = ecc_generate_key(ecc);

    mpz_init(mx);
    mpz_set_str(mx, "10", 10);
    mpz_init(my);
    mpz_set_str(my, "20", 10);
    Point message = point_new(mx, my);

    gmp_printf("public_key = (%Zd, %Zd)\n", key.public_key.x, key.public_key.y);
    gmp_printf("private_key = %Zd\n", key.private_key);
    gmp_printf("message = (%Zd, %Zd)\n\n", message.x, message.y);

    Ciphertext ciphertext = ecc_encrypt(ecc, &message, &key.public_key);
    gmp_printf("(P1.x, P1.y) = (%Zd, %Zd)\n", ciphertext.P1.x, ciphertext.P1.y);
    gmp_printf("(P2.x, P2.y) = (%Zd, %Zd)\n\n", ciphertext.P2.x, ciphertext.P2.y);

    Point decrypted = ecc_decrypt(ecc, &ciphertext, key.private_key);
    gmp_printf("decrypted = (%Zd, %Zd)\n", decrypted.x, decrypted.y);

    point_free(&decrypted);
    point_free(&ciphertext.P1);
    point_free(&ciphertext.P2);
    key_free(&key);
    point_free(&message);
    point_free(&G);
    ecc_free(ecc);
    elliptic_curve_free(curve);

    mpz_clear(prime);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(mx);
    mpz_clear(my);

    return 0;
}