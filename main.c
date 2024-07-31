#include "F28x_Project.h"
#include "driverlib.h"
#include "device.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <flecc_in_c/gfp/gfp_const_runtime.h>
#include <flecc_in_c/gfp/gfp.h>
#include <flecc_in_c/io/io.h>
#include <flecc_in_c/bi/bi.h>
#include <flecc_in_c/utils/param.h>
#include <flecc_in_c/utils/parse.h>

// Curve P-256
const char *curve_type_str = "secp256r1";
// R squared
const char *Rsq_str = "00000004FFFFFFFDFFFFFFFFFFFFFFFEFFFFFFFBFFFFFFFF0000000000000003";
// Coordinates of P
const char *p0x_str = "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
const char *p0y_str = "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
const char *p0z_str = "0000000000000000000000000000000000000000000000000000000000000001";

// Coordinates of Q
const char *q0x_str = "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
const char *q0y_str = "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
const char *q0z_str = "0000000000000000000000000000000000000000000000000000000000000001";

// k scalar value
const char *KB =  "1001101101011111110111"; //26D7F7


eccp_parameters_t curve_params;
gfp_prime_data_t prime_data;

// Registers
gfp_t R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10;
gfp_t T0, T1, T2, T3, T4, T5, T6;

gfp_t Qx, Qy, Qz;
gfp_t Px, Py, Pz;

gfp_t r_sq;
gfp_t X_A, Y_A;

int counter = 0;

void parse_bigint( const char *string, uint_t *big_int, const int bi_length ) {
    int len = strlen( string );
    bigint_parse_hex_var( big_int, bi_length, string, len );
}

void multiply(gfp_t res, gfp_t a, gfp_t b, const gfp_prime_data_t *prime_data) {
    gfp_cr_mont_multiply_sos(res, a, b, prime_data); // a * b * R^(-1)
    gfp_cr_mont_multiply_sos(res, res, r_sq, prime_data); // [a * b * R^(-1)] * R^2 * R^(-1)
}

/**
 *  MNAMNAA-based point doubling. Result = 2Q.
 *  @param T1 the x-coordinate of Q
 *  @param T2 the y-coordinate of Q
 *  @param T3 the z-coordinate of Q
 *  @param prime_data elliptic curve parameters
 */
void point_doubling(gfp_t T1, gfp_t T2, gfp_t T3, const gfp_prime_data_t *prime_data){
    // Delta 1
    // OP1: T4 = T3 * T3
    multiply(T4, T3, T3, prime_data);
    // * N
    gfp_cr_negate(T0, T1, prime_data);
    // T5 = T1 + T4
    gfp_cr_add(T5, T1, T4, prime_data);
    // T6 = T2 * T2
    multiply(T6, T2, T2, prime_data);
    // T4 = -T4
    gfp_cr_negate(T4, T4, prime_data);
    // T2 = T2 + T2
    gfp_cr_add(T2, T2, T2, prime_data);
    // T4 = T1 + T4
    gfp_cr_add(T4, T1, T4, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 2
    // OP8: T5 = T4 * T5
    multiply(T5, T4, T5, prime_data);
    // * N
    gfp_cr_negate(T0, T2, prime_data);
    // T4 = T5 + T5
    gfp_cr_add(T4, T5, T5, prime_data);
    // T3 = T2 * T3
    multiply(T3, T2, T3, prime_data);
    // * N
    gfp_cr_negate(T0, T4, prime_data);
    // T4 = T4 + T5
    gfp_cr_add(T4, T4, T5, prime_data);
    // T2 = T6 + T6
    gfp_cr_add(T2, T6, T6, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 3
    // OP 15: T5 = T4 * T4
    multiply(T5, T4, T4, prime_data);
    // * N
    gfp_cr_negate(T0, T1, prime_data);
    // T6 = T2 + T2
    gfp_cr_add(T6, T2, T2, prime_data);
    // T6 = T1 * T6
    multiply(T6, T1, T6, prime_data);
    // T1 = -T6
    gfp_cr_negate(T1, T6, prime_data);
    // T1 = T1 + T1
    gfp_cr_add(T1, T1, T1, prime_data);
    // T1 = T1 + T5
    gfp_cr_add(T1, T1, T5, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 4
    // OP 22: T2 = T2 * T2
    multiply(T2, T2, T2, prime_data);
    // T5 = -T1
    gfp_cr_negate(T5, T1, prime_data);
    // T5 = T5 + T6
    gfp_cr_add(T5, T5, T6, prime_data);
    // T5 = T4 * T5
    multiply(T5, T4, T5, prime_data);
    // T2 = -T2
    gfp_cr_negate(T2, T2, prime_data);
    // T2 = T2 + T2
    gfp_cr_add(T2, T2, T2, prime_data);
    // T2 = T2 + T5
    gfp_cr_add(T2, T2, T5, prime_data);

}


/**
 *  MNAMNAA-based point addition. Result = Q + P.
 *  @param X1 the x-coordinate of Q
 *  @param Y1 the y-coordinate of Q
 *  @param Z1 the z-coordinate of Q
 *  @param X the x-coordinate of P
 *  @param Y the y-coordinate of P
 *  @param prime_data elliptic curve parameters
 */
void point_addition(gfp_t X1, gfp_t Y1, gfp_t Z1, gfp_t X, gfp_t Y, const gfp_prime_data_t *prime_data){

    // Delta 1
    // OP 1: T4 = T3^2
    multiply(R4, Z1, Z1, prime_data);
    // * N
    gfp_cr_negate(R0, R1, prime_data);
    // * A
    gfp_cr_add(R5, R1, R4, prime_data);
    // T5 = Tx * T4
    multiply(R5, X, R4, prime_data);
    // T6 = -T1
    gfp_cr_negate(R6, X1, prime_data);
    // T5 = T5 + T6
    gfp_cr_add(R5, R5, R6, prime_data);
    // * A
    gfp_cr_add(R0, R1, R4, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 2
    // OP 8: T6 = T5^2
    multiply(R6, R5, R5, prime_data);
    // * N
    gfp_cr_negate(R0, R2, prime_data);
    // * A
    gfp_cr_add(R0, R5, R5, prime_data);
    // T7 = T1 * T6
    multiply(R7, X1, R6, prime_data);
    // * N
    gfp_cr_negate(R0, R4, prime_data);
    // T8 = T7 + T7
    gfp_cr_add(R8, R7, R7, prime_data);
    // * A
    gfp_cr_add(R0, R0, R0, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 3
    // OP 15: T9 = T5 * T6
    multiply(R9, R5, R6, prime_data);
    // * N
    gfp_cr_negate(R0, R1, prime_data);
    // T8 = T8 + T9
    gfp_cr_add(R8, R8, R9, prime_data);
    // T4 = T3 * T4
    multiply(R4, Z1, R4, prime_data);
    // * N
    gfp_cr_negate(R0, R6, prime_data);
    // * A
    gfp_cr_add(R0, R1, R1, prime_data);
    // * A
    gfp_cr_add(R0, R1, R5, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 4
    // OP 22: T4 = Ty * T4
    multiply(R4, Y, R4, prime_data);
    // T10 = -T2
    gfp_cr_negate(R10, Y1, prime_data);
    // T4 = T4 + T10
    gfp_cr_add(R4, R4, R10, prime_data);
    // T10 = T4 * T4
    multiply(R10, R4, R4, prime_data);
    // T8 = -T8
    gfp_cr_negate(R8, R8, prime_data);
    // T1 = T10 + T8
    gfp_cr_add(X1, R10, R8, prime_data);
    // * A
    gfp_cr_add(R0, R2, R5, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 5
    // OP 29: T8 = T2 * T9
    multiply(R8, Y1, R9, prime_data);
    // T6 = -T1
    gfp_cr_negate(R6, X1, prime_data);
    // T6 = T6 + T7
    gfp_cr_add(R6, R6, R7, prime_data);
    // T6 = T6 * T4
    multiply(R6, R6, R4, prime_data);
    // T9 = -T8
    gfp_cr_negate(R9, R8, prime_data);
    // T2 = T6 + T9
    gfp_cr_add(Y1, R6, R9, prime_data);
    // * A
    gfp_cr_add(R0, R2, R5, prime_data);

    // 2000 NOPs
    for( int k = 0; k < 2000; k++){
        counter++;
    }

    // Delta 6
    // OP 36: T3 = T3 * T5
    multiply(Z1, Z1, R5, prime_data);
    // T4 = -T7
    gfp_cr_negate(R4, R7, prime_data);
    // T4 = T1 + T4
    gfp_cr_add(R4, R1, R4, prime_data);
    // T5 = T4 * T4
    multiply(R5, R4, R4, prime_data);
    // T6 = -T8
    gfp_cr_negate(R6, R8, prime_data);
    // T6 = T2 + T6
    gfp_cr_add(R6, R2, R6, prime_data);
    // * A
    gfp_cr_add(R0, R2, R5, prime_data);

}


int main(void) {

    // Load EC P-256 curve
    curve_params.curve_type = param_get_curve_type_from_name( curve_type_str, strlen( curve_type_str ) );
    param_load( &curve_params, curve_params.curve_type );

    prime_data = curve_params.prime_data;
    //prime_data = curve_params.order_n_data;
    const gfp_prime_data_t * prime = &prime_data;
    //dataWords = curve_params.order_n_data.words;
    parse_bigint( Rsq_str, r_sq, curve_params.order_n_data.words);

    // Step 1 - Initially assign Q = P
    parse_bigint( q0x_str, Qx, curve_params.order_n_data.words); // assign q0x_str as value to Qx
    parse_bigint( q0y_str, Qy, curve_params.order_n_data.words);
    parse_bigint( q0z_str, Qz, curve_params.order_n_data.words);

    // P should not change its value
    parse_bigint( p0x_str, Px, curve_params.order_n_data.words); // assign p0x_str as value to Px
    parse_bigint( p0y_str, Py, curve_params.order_n_data.words);
    parse_bigint( p0z_str, Pz, curve_params.order_n_data.words);

    // Step2
    int l = strlen(KB);
    int t = l - 1;
    int i = l;

    for ( i = t; i > 0; i--) {
        // 5000 NOPs
        for( int k = 0; k < 5000; k++){
            counter++;
        }
        //Step3
        point_doubling(Qx, Qy, Qz, prime);

        // 5000 NOPs
        for( int j = 0; j < 5000; j++){
            counter++;
        }

        if (KB[l - i] == '1') {
            point_addition(Qx, Qy, Qz, Px, Py, prime);

        }
    }

    // Convert Jacobian coordinates to Affine coordinates
    if (KB[l - 1] == '1') {
        // X3^2
        multiply(R0, Qz, Qz, prime);
        bigint_copy_var(R1, R0, curve_params.order_n_data.words );
        // X_A = 1/X3^2
        gfp_binary_euclidean_inverse( X_A, R1, prime);
        // X_A = X1/X3^2
        multiply(X_A, Qx, X_A, prime);
        // X3^3
        multiply(R2, Qz, R0, prime);
        // Y_A = 1/X3^3
        gfp_binary_euclidean_inverse( Y_A, R2, prime );
        // Y_A = X2/X3^2
        multiply(Y_A, Qy, Y_A, prime);

        printf("\nX_A: ");
        io_print_bigint_var(X_A, curve_params.order_n_data.words);
        printf("Y_A: ");
        io_print_bigint_var(Y_A, curve_params.order_n_data.words);

    } else {
        printf("Error: The first bit of the key must be '1' for this algorithm to work!");
    }

    return 0;
}

/* Expected results */
// Using 26D7F7 as key:
// X_A: FDC45A02 2AB30C31 AE94ACEA 99AF39D8 D37E83A8 89C2EEE9 B5603DC9 9FD7B1DE
// Y_A: BF76B9A8 4B15F230 A663DF25 238DA790 156DE1C6 B72F17E7 7B1186BE AC840F40
