#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <openssl/sha.h>
#include <openssl/rand.h>
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

void ecbsa_generate_sbox(const Point* S, uint8_t sbox[256]) {
    for (int i = 0; i < 256; i++) sbox[i] = i;

    unsigned char hash[SHA256_DIGEST_LENGTH];
    char* x_str = mpz_get_str(NULL, 10, S->x);
    SHA256((unsigned char*)x_str, strlen(x_str), hash);
    free(x_str);

    unsigned char buf[4];
    for (int i = 255; i > 0; i--) {
        RAND_bytes(buf, sizeof(buf));
        int j = ((unsigned int)(buf[0] ^ hash[i % SHA256_DIGEST_LENGTH]) << 24 |
                 (buf[1] ^ hash[(i + 1) % SHA256_DIGEST_LENGTH]) << 16 |
                 (buf[2] ^ hash[(i + 2) % SHA256_DIGEST_LENGTH]) << 8 |
                 (buf[3] ^ hash[(i + 3) % SHA256_DIGEST_LENGTH])) % (i + 1);
        uint8_t temp = sbox[i];
        sbox[i] = sbox[j];
        sbox[j] = temp;
    }
}


void ecbsa_generate_inv_sbox(const uint8_t sbox[256], uint8_t inv_sbox[256]) {
    for (int i = 0; i < 256; i++) inv_sbox[i] = 0;
    for (int i = 0; i < 256; i++) inv_sbox[sbox[i]] = (uint8_t)i;
}

void ecbsa_key_schedule(ECC* ecc, const mpz_t d, uint8_t round_keys[ECBSA_ROUNDS + 1][ECBSA_BLOCK_SIZE]) {
    mpz_t dr;
    mpz_init(dr);
    for (int r = 0; r <= ECBSA_ROUNDS; ++r) {
        mpz_add_ui(dr, d, r);
        mpz_mod(dr, dr, ecc->curve->p);
        Point R = elliptic_curve_scalar_multiply(ecc->curve, &ecc->G, dr);
        mpz_t tmp;
        mpz_init_set(tmp, R.x);
        for (int i = 0; i < ECBSA_BLOCK_SIZE; i++) {
            round_keys[r][ECBSA_BLOCK_SIZE - 1 - i] = (uint8_t) mpz_fdiv_ui(tmp, 256);
            mpz_fdiv_q_ui(tmp, tmp, 256);
        }
        point_free(&R);
        mpz_clear(tmp);
    }
    mpz_clear(dr);
}

void shift_rows(uint8_t state[ECBSA_BLOCK_SIZE]) {
    uint8_t temp[ECBSA_BLOCK_SIZE];
    memcpy(temp, state, ECBSA_BLOCK_SIZE);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            state[i * 4 + j] = temp[i * 4 + ((j + i) % 4)];
        }
    }
}

static uint8_t gf_mult(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (int i = 0; i < 8; i++) {
        if (b & 1) p ^= a;
        uint8_t hi_bit = a & 0x80;
        a <<= 1;
        if (hi_bit) a ^= 0x1B;
        b >>= 1;
    }
    return p;
}

void mix_columns(uint8_t state[ECBSA_BLOCK_SIZE]) {
    uint8_t temp[4];
    for (int i = 0; i < 4; i++) {
        temp[0] = gf_mult(state[i * 4], 2) ^ gf_mult(state[i * 4 + 1], 3) ^ state[i * 4 + 2] ^ state[i * 4 + 3];
        temp[1] = state[i * 4] ^ gf_mult(state[i * 4 + 1], 2) ^ gf_mult(state[i * 4 + 2], 3) ^ state[i * 4 + 3];
        temp[2] = state[i * 4] ^ state[i * 4 + 1] ^ gf_mult(state[i * 4 + 2], 2) ^ gf_mult(state[i * 4 + 3], 3);
        temp[3] = gf_mult(state[i * 4], 3) ^ state[i * 4 + 1] ^ state[i * 4 + 2] ^ gf_mult(state[i * 4 + 3], 2);
        for (int j = 0; j < 4; j++) state[i * 4 + j] = temp[j];
    }
}

void inv_mix_columns(uint8_t state[ECBSA_BLOCK_SIZE]) {
    uint8_t temp[4];
    for (int i = 0; i < 4; i++) {
        temp[0] = gf_mult(state[i * 4], 0x0E) ^ gf_mult(state[i * 4 + 1], 0x0B) ^ gf_mult(state[i * 4 + 2], 0x0D) ^ gf_mult(state[i * 4 + 3], 0x09);
        temp[1] = gf_mult(state[i * 4], 0x09) ^ gf_mult(state[i * 4 + 1], 0x0E) ^ gf_mult(state[i * 4 + 2], 0x0B) ^ gf_mult(state[i * 4 + 3], 0x0D);
        temp[2] = gf_mult(state[i * 4], 0x0D) ^ gf_mult(state[i * 4 + 1], 0x09) ^ gf_mult(state[i * 4 + 2], 0x0E) ^ gf_mult(state[i * 4 + 3], 0x0B);
        temp[3] = gf_mult(state[i * 4], 0x0B) ^ gf_mult(state[i * 4 + 1], 0x0D) ^ gf_mult(state[i * 4 + 2], 0x09) ^ gf_mult(state[i * 4 + 3], 0x0E);
        for (int j = 0; j < 4; j++) state[i * 4 + j] = temp[j];
    }
}

void inv_shift_rows(uint8_t state[ECBSA_BLOCK_SIZE]) {
    uint8_t temp[ECBSA_BLOCK_SIZE];
    memcpy(temp, state, ECBSA_BLOCK_SIZE);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            state[i * 4 + j] = temp[i * 4 + ((j - i + 4) % 4)];
        }
    }
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
    uint8_t round_keys[ECBSA_ROUNDS + 1][ECBSA_BLOCK_SIZE];
    ecbsa_key_schedule(ecc, d, round_keys);

    uint8_t state[ECBSA_BLOCK_SIZE];
    memcpy(state, in, ECBSA_BLOCK_SIZE);

    for (int r = 0; r < ECBSA_ROUNDS; r++) {
        for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
            state[i] = sbox[state[i]];
        shift_rows(state);
        if (r < ECBSA_ROUNDS - 1)
            mix_columns(state);
        for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
            state[i] ^= round_keys[r][i];
    }
    for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
        state[i] = sbox[state[i]];
    shift_rows(state);
    for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
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
    uint8_t round_keys[ECBSA_ROUNDS + 1][ECBSA_BLOCK_SIZE];
    ecbsa_generate_inv_sbox(sbox, inv_sbox);
    ecbsa_key_schedule(ecc, d, round_keys);

    uint8_t state[ECBSA_BLOCK_SIZE];
    memcpy(state, in, ECBSA_BLOCK_SIZE);

    for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
        state[i] ^= round_keys[ECBSA_ROUNDS][i];
    inv_shift_rows(state);
    for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
        state[i] = inv_sbox[state[i]];

    for (int r = ECBSA_ROUNDS - 1; r >= 0; r--) {
        for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
            state[i] ^= round_keys[r][i];
        if (r < ECBSA_ROUNDS - 1)
            inv_mix_columns(state);
        inv_shift_rows(state);
        for (int i = 0; i < ECBSA_BLOCK_SIZE; i++)
            state[i] = inv_sbox[state[i]];
    }

    memcpy(out, state, ECBSA_BLOCK_SIZE);
}

void compute_d_from_Sx(const Point* S, const mpz_t p, mpz_t d) {
    unsigned char hash[SHA256_DIGEST_LENGTH];
    char* x_str = mpz_get_str(NULL, 10, S->x);
    SHA256((unsigned char*)x_str, strlen(x_str), hash);
    free(x_str);
    mpz_import(d, SHA256_DIGEST_LENGTH, 1, 1, 0, 0, hash);
    mpz_mod(d, d, p);
}