#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "elliptic_curve.h"
#include "ecc.h"
#include "ecbsa.h"

int main() {
    mpz_t p, a, b, Gx, Gy, d;
    mpz_init_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_init_set_str(a, "0", 10);
    mpz_init_set_str(b, "7", 10);
    mpz_init_set_str(Gx, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_init_set_str(Gy, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_init_set_str(d, "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);

    Elliptic_Curve* curve = elliptic_curve_new(a,b,p);
    Point G = point_new(Gx,Gy);
    mpz_mod(G.x, G.x, curve->p);
    mpz_mod(G.y, G.y, curve->p);
    ECC* ecc = ecc_new(curve, &G);

    const char* msg = "Hello, World!";
    size_t msg_len = strlen(msg);
    size_t block_size = ECBSA_BLOCK_SIZE;
    size_t padded_len = ((msg_len + block_size - 1) / block_size) * block_size;

    uint8_t* plaintext = malloc(padded_len);
    memcpy(plaintext, msg, msg_len);
    add_padding(plaintext, msg_len, block_size);

    uint8_t* ciphertext = malloc(padded_len);
    uint8_t* recovered = malloc(padded_len);
    uint8_t sbox[256];

    ecbsa_generate_sbox(ecc, sbox);

    ecbsa_encrypt(ecc, d, sbox, plaintext, msg_len, ciphertext);

    printf("Ciphertext: ");
    for (size_t i = 0; i < padded_len; ++i) printf("%02X ", ciphertext[i]);
    printf("\n");

    size_t recovered_len = 0;
    ecbsa_decrypt(ecc, d, sbox, ciphertext, padded_len, recovered, &recovered_len);

    printf("Recovered: ");
    for (size_t i = 0; i < recovered_len; ++i) printf("%02X ", recovered[i]);
    printf("\n");

    char* outmsg = malloc(recovered_len + 1);
    memcpy(outmsg, recovered, recovered_len);
    outmsg[recovered_len] = '\0';
    printf("Recovered string: %s\n", outmsg);

    free(plaintext);
    free(ciphertext);
    free(recovered);
    free(outmsg);
    mpz_clears(p,a,b,Gx,Gy,d, NULL);
    point_free(&G);
    ecc_free(ecc);
    elliptic_curve_free(curve);
    return 0;
}