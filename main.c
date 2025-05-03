#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "elliptic_curve.h"
#include "ecc.h"
#include "ecbsa.h"

int main() {
    // Initialize elliptic curve parameters
    mpz_t p, a, b, Gx, Gy;
    mpz_init_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_init_set_str(a, "0", 10);
    mpz_init_set_str(b, "7", 10);
    mpz_init_set_str(Gx, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_init_set_str(Gy, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

    Elliptic_Curve* curve = elliptic_curve_new(a, b, p);
    Point G = point_new(Gx, Gy);
    mpz_mod(G.x, G.x, curve->p);
    mpz_mod(G.y, G.y, curve->p);
    ECC* ecc = ecc_new(curve, &G);

    // Sender: Generate private key k_A
    mpz_t k_A;
    mpz_init(k_A);
    generate_random_private_key(k_A, p);
    printf("Sender's private key k_A: ");
    mpz_out_str(stdout, 16, k_A);
    printf("\n");

    // Sender: Generate public key P_A = [k_A]G
    Point P_A = elliptic_curve_scalar_multiply(curve, &G, k_A);
    printf("Sender's public key P_A: (");
    mpz_out_str(stdout, 16, P_A.x);
    printf(", ");
    mpz_out_str(stdout, 16, P_A.y);
    printf(")\n\n");

    // Receiver: Generate private key k_B
    mpz_t k_B;
    mpz_init(k_B);
    generate_random_private_key(k_B, p);
    printf("Receiver's private key k_B: ");
    mpz_out_str(stdout, 16, k_B);
    printf("\n");

    // Receiver: Generate public key P_B = [k_B]G
    Point P_B = elliptic_curve_scalar_multiply(curve, &G, k_B);
    printf("Receiver's public key P_B: (");
    mpz_out_str(stdout, 16, P_B.x);
    printf(", ");
    mpz_out_str(stdout, 16, P_B.y);
    printf(")\n\n");

    // Sender: Compute shared secret S_sender = [k_A]P_B
    Point S_sender = elliptic_curve_scalar_multiply(curve, &P_B, k_A);
    printf("Sender's shared secret S_sender: (");
    mpz_out_str(stdout, 16, S_sender.x);
    printf(", ");
    mpz_out_str(stdout, 16, S_sender.y);
    printf(")\n");

    // Receiver: Compute shared secret S_receiver = [k_B]P_A
    Point S_receiver = elliptic_curve_scalar_multiply(curve, &P_A, k_B);
    printf("Receiver's shared secret S_receiver: (");
    mpz_out_str(stdout, 16, S_receiver.x);
    printf(", ");
    mpz_out_str(stdout, 16, S_receiver.y);
    printf(")\n\n");

    // Sender: Compute d_sender = SHA-256(S_sender.x) mod p
    mpz_t d_sender;
    mpz_init(d_sender);
    compute_d_from_Sx(&S_sender, p, d_sender);
    printf("Sender's d_sender: ");
    mpz_out_str(stdout, 16, d_sender);
    printf("\n");

    // Receiver: Compute d_receiver = SHA-256(S_receiver.x) mod p
    mpz_t d_receiver;
    mpz_init(d_receiver);
    compute_d_from_Sx(&S_receiver, p, d_receiver);
    printf("Receiver's d_receiver: ");
    mpz_out_str(stdout, 16, d_receiver);
    printf("\n\n");

    // Generate S-box
    uint8_t sbox[256];
    ecbsa_generate_sbox(&S_receiver, sbox);
    printf("S-box generated\n\n");

    // Encryption (Sender)
    const char* msg = "Hello, World!";
    size_t msg_len = strlen(msg);
    size_t block_size = ECBSA_BLOCK_SIZE;
    size_t padded_len = ((msg_len + block_size - 1) / block_size) * block_size;
    uint8_t* plaintext = malloc(padded_len);
    memcpy(plaintext, msg, msg_len);
    add_padding(plaintext, msg_len, block_size);
    uint8_t* ciphertext = malloc(padded_len);
    ecbsa_encrypt(ecc, d_sender, sbox, plaintext, msg_len, ciphertext);
    printf("Encrypted message: ");
    for (size_t i = 0; i < padded_len; ++i) printf("%02X ", ciphertext[i]);
    printf("\n\n");

    // Decryption (Receiver)
    uint8_t* recovered = malloc(padded_len);
    size_t recovered_len = 0;
    ecbsa_decrypt(ecc, d_receiver, sbox, ciphertext, padded_len, recovered, &recovered_len);
    printf("Decrypted message (before inverse S-box): ");
    for (size_t i = 0; i < padded_len; ++i) printf("%02X ", recovered[i]);
    printf("\n\n");

    // Output final decrypted message
    printf("Decrypted message: %.*s\n", (int)recovered_len, recovered);

    // Cleanup
    free(plaintext);
    free(ciphertext);
    free(recovered);
    mpz_clears(p, a, b, Gx, Gy, k_A, k_B, d_sender, d_receiver, NULL);
    point_free(&G);
    point_free(&P_A);
    point_free(&P_B);
    point_free(&S_sender);
    point_free(&S_receiver);
    ecc_free(ecc);
    elliptic_curve_free(curve);
    return 0;
}