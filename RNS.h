/*
 ============================================================================
 Name        : RNS_ECC.c
 Author      : Apostolos P. Fournaris
 Version     :
 Copyright   : Your copyright notice
 Description : SCA-FA RNS Scalar multiplier
 ============================================================================
 */

#ifndef RNS_H_
#define RNS_H_

#define MOD_NUM 4
#define RAND_COMB_NUM 70  //combinations of MOD_NUM out of 2*MOD_NUM
#define EC_Field_P_bit_lenght 192

#include "gmp.h"
#include <sys/time.h>

typedef struct {
	mpz_t m0;
	mpz_t m1;
	mpz_t m2;
	mpz_t m3;
}rns_num;

typedef struct {
	mpz_t MM;
	mpz_t m[2*MOD_NUM];
	}MmodP;

typedef struct {
	int m0_num;
	int m1_num;
	int m2_num;
	int m3_num;
	int moduli[MOD_NUM];
	int base_index;
	mpz_t Mb;
	int m_nn_nums[MOD_NUM];
}base;



typedef struct {
	mpz_t modulo;
	mpz_t p_mod;
	mpz_t min_inv_p_mod;
	mpz_t inv_modulo;
	mpz_t rel_inv_mod_matr[2*MOD_NUM];
	mpz_t Mb_inv_mi[RAND_COMB_NUM];
}rns_base_element_data[2*MOD_NUM];

typedef struct {
	mpz_t x[2*MOD_NUM];
	mpz_t y[2*MOD_NUM];
	mpz_t z[2*MOD_NUM];
	}ec_point_rns;

typedef struct {
	mpz_t x;
	mpz_t y;
	mpz_t z;
}ec_point;

int Bn[MOD_NUM+1];
int Bnn[MOD_NUM+1];

rns_base_element_data B_pool;

base base_index[RAND_COMB_NUM];

ec_point_rns rns_EC_point_init(void);
rns_num rns_num_init(void);
ec_point EC_point_init(void);
//rns_base_element_data rns_base_data_init(void);
void generate_base_element_pool(rns_base_element_data pool, mpz_t field_p, base base_rand_index[],MmodP *Mp);
void Base_extention(rns_num num_b1, rns_num num_b2, int base_1[], int base_2[], rns_base_element_data pool);
void binary_to_rns(mpz_t num, mpz_t num_b[], int base[], rns_base_element_data pool);
void rns_to_binary(mpz_t *num, mpz_t num_b1[], int base[], rns_base_element_data pool);
void Base_extention_B(mpz_t num_b1[], mpz_t num_b2[], int base_1[], int base_2[], rns_base_element_data pool);
void RNS_mul(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool);
void RNS_add(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool);
void RNS_sub(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool);
void RNS_Add_mod_M(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool);
void RNS_Sub_mod(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool, mpz_t field_p);
void RNS_Add_mod(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool, mpz_t field_p);
void RNS_add_no(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool);
void RNS_sub_no(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool);
void RNS_Montg_mul(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool);
void EdC_point_add(ec_point_rns *P3, ec_point_rns P1, ec_point_rns P2, int base1[], int base2[], rns_base_element_data pool, mpz_t field_p);
void EdC_point_double(ec_point_rns *P3, ec_point_rns P1, int base1[], int base2[], rns_base_element_data pool, mpz_t field_p);
void EC_scalar_mul(ec_point *Q, ec_point P, ec_point R, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p);
void EC_scalar_mul_no_Base_Perm(ec_point *Q, ec_point P, ec_point R, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p);
void EC_scalar_mul_no_R_point(ec_point *Q, ec_point P, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p);
void EC_scalar_mul_no_Base_Perm_no_R_point(ec_point *Q, ec_point P, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p);
void Rand_EC_point(ec_point *Q, ec_point base_pointP, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p);
void rns_num_clear(rns_num t);
void ec_point_clear(ec_point t);
void ec_point_rns_clear(ec_point_rns t);

#endif /* RNS_H_ */
