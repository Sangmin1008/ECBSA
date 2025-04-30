#ifndef ECBSA_H
#define ECBSA_H

#include <stdint.h>
#include <gmp.h>
#include "ecc.h"

#define ECBSA_BLOCK_SIZE 16
#define ECBSA_ROUNDS 4

void add_padding(uint8_t* data, size_t len, size_t block_size);
void remove_padding(uint8_t* data, size_t* len);

void ecbsa_generate_sbox(ECC* ecc, uint8_t sbox[256]);
void ecbsa_generate_inv_sbox(const uint8_t sbox[256], uint8_t inv_sbox[256]);
void ecbsa_key_schedule(ECC* ecc, const mpz_t d, uint8_t round_keys[ECBSA_ROUNDS+1][ECBSA_BLOCK_SIZE]);

void ecbsa_encrypt_block(ECC* ecc, const mpz_t d, const uint8_t sbox[256], const uint8_t in[ECBSA_BLOCK_SIZE], uint8_t out[ECBSA_BLOCK_SIZE]);
void ecbsa_decrypt_block(ECC* ecc, const mpz_t d, const uint8_t sbox[256], const uint8_t in[ECBSA_BLOCK_SIZE], uint8_t out[ECBSA_BLOCK_SIZE]);

void ecbsa_encrypt(ECC* ecc, const mpz_t d, const uint8_t sbox[256],
                   const uint8_t* plaintext, size_t plaintext_len,
                   uint8_t* ciphertext);
void ecbsa_decrypt(ECC* ecc, const mpz_t d, const uint8_t sbox[256],
                   const uint8_t* ciphertext, size_t ciphertext_len,
                   uint8_t* recovered, size_t* recovered_len);

#endif // ECBSA_H