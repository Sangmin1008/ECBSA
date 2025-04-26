#include <stdio.h>
#include <time.h>
#include <C:\msys64\gmp-6.3.0\include\gmp.h>
#include <stdlib.h>
#include "elliptic_curve.h"
#include "ecc.h"

int main() {
    srand(time(NULL));
    mpz_t prime, a, b, x, y;
    mpz_init_set_ui(prime, 97);
    mpz_init_set_ui(a, 2);
    mpz_init_set_ui(b, 3);
    mpz_init_set_ui(x, 3);
    mpz_init_set_ui(y, 6);

    Elliptic_Curve* curve = elliptic_curve_new(a, b, prime);
    Point G = point_new(x, y);
    ECC* ecc = ecc_new(curve, &G);

    Key key = ecc_generate_key(ecc);
    mpz_init_set_ui(x, 10);
    mpz_init_set_ui(y, 20);
    Point message = point_new(x, y);

    gmp_printf("public_key = (%Zd, %Zd)\n", key.public_key.x, key.public_key.y);
    gmp_printf("private_key = %Zd\n", key.private_key);
    gmp_printf("message = (%Zd, %Zd)\n\n", message.x, message.y);

    Ciphertext ciphertext = ecc_encrypt(ecc, &message, &key.public_key);
    gmp_printf("(P1.x, P1.y) = (%Zd, %Zd)\n", ciphertext.P1.x, ciphertext.P1.y);
    gmp_printf("(P2.x, P2.y) = (%Zd, %Zd)\n\n", ciphertext.P2.x, ciphertext.P2.y);

    Point decrypted = ecc_decrypt(ecc, &ciphertext, key.private_key);
    gmp_printf("decrypted = (%Zd, %Zd)\n", decrypted.x, decrypted.y);

    point_free(&decrypted);
    point_free(&ciphertext.P1); point_free(&ciphertext.P2);
    key_free(&key);
    point_free(&message);
    point_free(&G);
    ecc_free(ecc);
    elliptic_curve_free(curve);
    mpz_clear(prime); mpz_clear(a); mpz_clear(b); mpz_clear(x); mpz_clear(y);

    return 0;
}

