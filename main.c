#include <stdio.h>
#include "elliptic_curve.h"
#include "ecc.h"

int main() {
    int64_t prime = 97;
    Elliptic_Curve* curve = elliptic_curve_new(2, 3, prime);
    Point G = point_new(3, 6);
    ECC* ecc = ecc_new(curve, &G);

    Key key = ecc_generate_key(ecc);
    Point message = point_new(10, 20);
    printf("public_key = %d\nprivate_key = %d\n", key.public_key, key.private_key);
    printf("message = (%d, %d)\n\n", message.x, message.y);

    Ciphertext ciphertext = ecc_encrypt(ecc, &message, &key.public_key);
    printf("(P1.x, P1.y) = (%d, %d)\n", ciphertext.P1.x, ciphertext.P1.y);
    printf("(P2.x, P2.y) = (%d, %d)\n\n", ciphertext.P2.x, ciphertext.P2.y);
    Point decrypted = ecc_decrypt(ecc, &ciphertext, key.private_key);

    printf("decrypted = (%d, %d)", decrypted.x, decrypted.y);
}