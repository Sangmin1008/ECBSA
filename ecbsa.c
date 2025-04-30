#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "ecbsa.h"

void add_padding(uint8_t* data, size_t len, size_t block_size) {
    uint8_t pad_value = block_size - (len % block_size);
    for (size_t i = len; i < len + pad_value; ++i) {
        data[i] = pad_value;
    }
}

void remove_padding(uint8_t* data, size_t* len) {
    uint8_t pad_value = data[*len - 1];
    *len -= pad_value;
}

void ecbsa_generate_sbox(ECC* ecc, uint8_t sbox[256]) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    mpz_t k;
    mpz_init(k);
    mpz_urandomm(k, state, ecc->curve->p);

    Point P = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, k);
    char* x_str = mpz_get_str(NULL, 10, P.x);
    unsigned long seed = 0;
    for (int i = 0; x_str[i] != '\0'; ++i) {
        seed = (seed * 31 + x_str[i]) % ULONG_MAX;
    }
    free(x_str);
    srand(seed);

    for (int i = 0; i < 256; ++i) sbox[i] = i;
    for (int i = 255; i > 0; --i) {
        int j = rand() % (i + 1);
        uint8_t temp = sbox[i];
        sbox[i] = sbox[j];
        sbox[j] = temp;
    }

    point_free(&P);
    mpz_clear(k);
    gmp_randclear(state);
}

void ecbsa_generate_inv_sbox(const uint8_t sbox[256], uint8_t inv_sbox[256]) {
    for (int i = 0; i < 256; ++i) inv_sbox[i] = 0;
    for (int i = 0; i < 256; ++i) inv_sbox[sbox[i]] = (uint8_t)i;
}

void ecbsa_key_schedule(ECC* ecc, const mpz_t d, uint8_t round_keys[ECBSA_ROUNDS+1][ECBSA_BLOCK_SIZE]) {
    mpz_t dr;
    mpz_init(dr);
    for (int r = 0; r <= ECBSA_ROUNDS; ++r) {
        mpz_add_ui(dr, d, r);
        mpz_mod(dr, dr, ecc->curve->p);
        Point R = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, dr);
        mpz_t tmp;
        mpz_init_set(tmp, R.x);
        for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i) {
            round_keys[r][ECBSA_BLOCK_SIZE-1-i] = (uint8_t) mpz_fdiv_ui(tmp, 256);
            mpz_fdiv_q_ui(tmp, tmp, 256);
        }
        point_free(&R);
        mpz_clear(tmp);
    }
    mpz_clear(dr);
}

void ecbsa_encrypt(ECC* ecc, const mpz_t d, const uint8_t sbox[256],
                   const uint8_t* plaintext, size_t plaintext_len,
                   uint8_t* ciphertext) {
    size_t block_size = 16;
    size_t padded_len = ((plaintext_len + block_size - 1) / block_size) * block_size;
    uint8_t* padded_plaintext = malloc(padded_len);

    memcpy(padded_plaintext, plaintext, plaintext_len);
    add_padding(padded_plaintext, plaintext_len, block_size);

    for (size_t i = 0; i < padded_len; i += block_size) {
        ecbsa_encrypt_block(ecc, d, sbox, padded_plaintext + i, ciphertext + i);
    }

    free(padded_plaintext);
}

void ecbsa_encrypt_block(ECC* ecc, const mpz_t d,
                         const uint8_t sbox[256],
                         const uint8_t in[ECBSA_BLOCK_SIZE],
                         uint8_t out[ECBSA_BLOCK_SIZE]) {
    uint8_t round_keys[ECBSA_ROUNDS+1][ECBSA_BLOCK_SIZE];
    ecbsa_key_schedule(ecc, d, round_keys);

    uint8_t state[ECBSA_BLOCK_SIZE];
    memcpy(state, in, ECBSA_BLOCK_SIZE);

    for (int r = 0; r < ECBSA_ROUNDS; ++r) {
        for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
            state[i] = sbox[state[i]];
        for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
            state[i] ^= round_keys[r][i];
    }
    for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
        state[i] = sbox[state[i]];
    for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
        state[i] ^= round_keys[ECBSA_ROUNDS][i];

    memcpy(out, state, ECBSA_BLOCK_SIZE);
}

void ecbsa_decrypt(ECC* ecc, const mpz_t d, const uint8_t sbox[256],
                   const uint8_t* ciphertext, size_t ciphertext_len,
                   uint8_t* recovered, size_t* recovered_len) {
    size_t block_size = 16;

    for (size_t i = 0; i < ciphertext_len; i += block_size) {
        ecbsa_decrypt_block(ecc, d, sbox, ciphertext + i, recovered + i);
    }

    *recovered_len = ciphertext_len;
    remove_padding(recovered, recovered_len);
}

void ecbsa_decrypt_block(ECC* ecc, const mpz_t d,
                         const uint8_t sbox[256],
                         const uint8_t in[ECBSA_BLOCK_SIZE],
                         uint8_t out[ECBSA_BLOCK_SIZE]) {
    uint8_t inv_sbox[256];
    uint8_t round_keys[ECBSA_ROUNDS+1][ECBSA_BLOCK_SIZE];
    ecbsa_generate_inv_sbox(sbox, inv_sbox);
    ecbsa_key_schedule(ecc, d, round_keys);

    uint8_t state[ECBSA_BLOCK_SIZE];
    memcpy(state, in, ECBSA_BLOCK_SIZE);

    for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
        state[i] ^= round_keys[ECBSA_ROUNDS][i];
    for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
        state[i] = inv_sbox[state[i]];

    for (int r = ECBSA_ROUNDS-1; r >= 0; --r) {
        for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
            state[i] ^= round_keys[r][i];
        for (int i = 0; i < ECBSA_BLOCK_SIZE; ++i)
            state[i] = inv_sbox[state[i]];
    }
    memcpy(out, state, ECBSA_BLOCK_SIZE);
}