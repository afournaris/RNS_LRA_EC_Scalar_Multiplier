/*
 ============================================================================
 Name        : RNS_ECC.c
 Author      : Apostolos P. Fournaris
 Version     :
 Copyright   : Your copyright notice
 Description : SCA-FA RSN Scalar multiplier
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include </home/user/workspaceA/RNS_ECC/RNS.h>

#include <gmp.h>
//#define DEBUG
//#define DEB_MONT



ec_point_rns rns_EC_point_init(void){
	ec_point_rns t;
	int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(t.x[i]);
		mpz_init(t.y[i]);
		mpz_init(t.z[i]);
	}
	return t;
}

/*rns_num rns_num_init(void){
	rns_num t;
	mpz_init(t.m0);
	mpz_init(t.m1);
	mpz_init(t.m2);
	mpz_init(t.m3);
	return t;
}*/

ec_point EC_point_init(void){
	ec_point t;
	mpz_init(t.x);
	mpz_init(t.y);
	mpz_init(t.z);
	return t;
}

/*rns_base_element_data rns_base_data_init(void){
	rns_base_element_data t;
		mpz_init(t.inv_modulo);
		mpz_init(t.min_inv_p_mod);
		mpz_init(t.p_mod);
		int i=0;
		for (i=0;i<2*MOD_NUM;i++){
			mpz_init(t.rel_inv_mod_matr[i]);
		}
		return t;
}*/

 void generate_base_element_pool(rns_base_element_data pool, mpz_t field_p, base base_rand_index[],MmodP *Mp){
	 int f;

	 for (f=0;f<2*MOD_NUM;f++){
			mpz_init(pool[f].inv_modulo);
			mpz_init(pool[f].min_inv_p_mod);
			mpz_init(pool[f].p_mod);
			int g=0;
			for (g=0;g<2*MOD_NUM;g++){
				mpz_init(pool[f].rel_inv_mod_matr[g]);
			}
		}
	 mpz_t tempMb;
	 mpz_init(tempMb);

//Here the 2 bases moduli are added. Currently supporting 4 moduli bases TODO autogenerate moduli
	 mpz_set_str(pool[0].modulo,"1125899905794047",10); // m0=2^50-2^20-1
	 mpz_set_str(pool[1].modulo,"1125899902648319",10); //m1=2^50-2^22-1
	 mpz_set_str(pool[2].modulo,"1125899906580479",10); //m2=2^50-2^18-1
	 mpz_set_str(pool[3].modulo,"1125899906841599",10); //m3=2^50-2^10-1
	 mpz_set_str(pool[4].modulo,"1125899906842624",10); //m4=2^50
	 mpz_set_str(pool[5].modulo,"1125899906842623",10);  //m5= 2^50-1
	 mpz_set_str(pool[6].modulo,"2251799813685247",10); //m6=2^51-1
	 mpz_set_str(pool[7].modulo,"562949953421311",10);	 //m7=2^49-1
//--------------------------------------------------------------------------------------------

	 int q;
		mpz_set(tempMb,pool[0].modulo);
	    for (q=1;q<2*MOD_NUM;q++){
			mpz_mul(tempMb, tempMb, pool[q].modulo);
	    }
		mpz_mod(tempMb,tempMb,field_p);
		for (q=0;q<2*MOD_NUM;q++){
			mpz_mod(Mp->m[q], tempMb, pool[q].modulo);
		}



// TODO The following code must be updated for autogeneration
	 int tt=0;
		 int i;
	 	 for (i=0; i<5;i++){

	 		 int j;

	 		 for (j=i+1;j<6;j++){

	 			 int k;

	 			 for (k=j+1;k<7;k++){
	 				 int l;

	 				 for (l=k+1;l<8;l++){
	 			 		base_rand_index[tt].base_index=tt;
	 			 		base_rand_index[tt].moduli[0]=i;
	 			 		base_rand_index[tt].moduli[1]=j;
						base_rand_index[tt].moduli[2]=k;
						base_rand_index[tt].moduli[3]=l;
	 			 		base_rand_index[tt].m0_num=i;
						base_rand_index[tt].m1_num=j;
	 			 		base_rand_index[tt].m2_num=k;
	 					base_rand_index[tt].m3_num=l;
	 					int t,s;
	 					s=0;
	 					for (t=0;t<8;t++){
	 						if (t!=i & t!=j & t!=k & t!=l) {base_rand_index[tt].m_nn_nums[s]=t; s=s+1;}
	 					}
	 					mpz_set(tempMb,pool[i].modulo);
	 					mpz_mul(tempMb, tempMb, pool[j].modulo);
	 					mpz_mul(tempMb, tempMb, pool[k].modulo);
	 					mpz_mul(tempMb, tempMb, pool[l].modulo);
	 					mpz_init_set(base_rand_index[tt].Mb,tempMb);

	 					 tt=tt+1;

	 				 }
	 			 }
	 		 }
	 	 }

//------------------------------------------------

	 mpz_t temp1,temp2,temp3;
	 mpz_init(temp1);
	 mpz_init(temp2);
	 mpz_init(temp3);



	 for (i=0;i<2*MOD_NUM;i++){
		 mpz_mod(pool[i].p_mod,field_p, pool[i].modulo);
		 mpz_neg(temp3,field_p);
		mpz_mod(temp3,temp3,pool[i].modulo);
		 mpz_gcdext (temp1, pool[i].min_inv_p_mod, temp2, temp3, pool[i].modulo);
		 mpz_mod(pool[i].min_inv_p_mod, pool[i].min_inv_p_mod, pool[i].modulo);
		 int j;
		 for (j=0;j<2*MOD_NUM;j++){
			 mpz_gcdext (temp1, pool[i].rel_inv_mod_matr[j], temp2, pool[i].modulo, pool[j].modulo);
		 }
		 j=0;
		 for(j=0;j<RAND_COMB_NUM;j++){
			 if (mpz_cmp(pool[base_rand_index[j].m0_num].modulo, pool[i].modulo)==0 || mpz_cmp(pool[base_rand_index[j].m1_num].modulo,pool[i].modulo)==0 || mpz_cmp(pool[base_rand_index[j].m2_num].modulo, pool[i].modulo)==0 ||mpz_cmp(pool[base_rand_index[j].m3_num].modulo,pool[i].modulo)==0){
				 mpz_init(pool[i].Mb_inv_mi[j]);
				 mpz_set_str(pool[i].Mb_inv_mi[j],"0",16);

				 }

			 else{
			 mpz_init(pool[i].Mb_inv_mi[j]);
			 mpz_gcdext (temp1, pool[i].Mb_inv_mi[j], temp2, base_rand_index[j].Mb, pool[i].modulo);
			 mpz_mod(pool[i].Mb_inv_mi[j], pool[i].Mb_inv_mi[j], pool[i].modulo);

		  }
		 }
	 }

	 mpz_clear(temp1);
	 mpz_clear(temp2);
	 mpz_clear(temp3);
	 mpz_clear(tempMb);

}




void binary_to_rns(mpz_t num, mpz_t num_b[], int base[], rns_base_element_data pool){
	int i=0;
	for (i=0;i<MOD_NUM;i++){
		mpz_mod(num_b[base[i]], num, pool[base[i]].modulo);
	}
}

void rns_to_binary(mpz_t *num, mpz_t num_b1[], int base[], rns_base_element_data pool){
mpz_t mrs_num[MOD_NUM];
int i=0;
int j=0;
for (i=0;i<MOD_NUM;i++){
	mpz_init(mrs_num[i]);
}

		//RNS to MRS conversion

		 mpz_set(mrs_num[0], num_b1[base[0]]);
		 for (i=1;i<MOD_NUM;i++){
			 mpz_sub(mrs_num[i],num_b1[base[i]],mrs_num[0]);
			 mpz_mul(mrs_num[i],mrs_num[i],pool[base[0]].rel_inv_mod_matr[base[i]]);
			 if (i==1){mpz_mod(mrs_num[i],mrs_num[i], pool[base[i]].modulo);}
			 else{
				 for (j=1;j<i;j++){
					 mpz_sub(mrs_num[i],mrs_num[i],mrs_num[j]);
					 mpz_mul(mrs_num[i],mrs_num[i],pool[base[j]].rel_inv_mod_matr[base[i]]);
					 mpz_mod(mrs_num[i],mrs_num[i], pool[base[i]].modulo);
				 }
			 }
		 }


		 //now mrs_num is in MRS format
		 //..............................................

		 //Conversion from MRS to binary
		 mpz_set(*num, mrs_num[MOD_NUM-1]);
		 for (i=MOD_NUM-2;i>=0;i--){
			 mpz_mul(*num,*num, pool[base[i]].modulo);
			 mpz_add(*num,*num,mrs_num[i]);
		 }
		 for (i=0;i<MOD_NUM;i++){
		 	mpz_clear(mrs_num[i]);
		 }

}

void Base_extention_B(mpz_t num_b1[], mpz_t num_b2[], int base_1[], int base_2[], rns_base_element_data pool){
	 mpz_t temp1;
	 mpz_init(temp1);
	 rns_to_binary(&temp1,num_b1,base_1,pool);
	 binary_to_rns(temp1, num_b2, base_2, pool);
	 mpz_clear(temp1);
}

void RNS_mul(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool){
	int i=0;
	for(i=0;i<MOD_NUM;i++){
		mpz_mul(num_r[base[i]],num_o1[base[i]],num_o2[base[i]]);
		mpz_mod(num_r[base[i]], num_r[base[i]],pool[base[i]].modulo);
	}
}

void RNS_add(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool){
	int i=0;
	for(i=0;i<MOD_NUM;i++){
		mpz_add(num_r[base[i]],num_o1[base[i]],num_o2[base[i]]);

		mpz_sub(num_r[base[i]],num_r[base[i]],pool[base[i]].p_mod);
		mpz_mod(num_r[base[i]], num_r[base[i]],pool[base[i]].modulo);
		}
}

void RNS_sub(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool ){
	int i=0;
	for(i=0;i<MOD_NUM;i++){
		mpz_sub(num_r[base[i]],num_o1[base[i]],num_o2[base[i]]);

		mpz_add(num_r[base[i]],num_r[base[i]],pool[base[i]].p_mod);
		mpz_mod(num_r[base[i]], num_r[base[i]],pool[base[i]].modulo);


	}
}

	void RNS_add_no(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool){
		int i=0;
		for(i=0;i<MOD_NUM;i++){
			mpz_add(num_r[base[i]],num_o1[base[i]],num_o2[base[i]]);
			mpz_mod(num_r[base[i]], num_r[base[i]],pool[base[i]].modulo);

		}


	}

	void RNS_sub_no(mpz_t num_r[], mpz_t num_o1[],mpz_t num_o2[], int base[], rns_base_element_data pool ){
		int i=0;
		for(i=0;i<MOD_NUM;i++){
			mpz_sub(num_r[base[i]],num_o1[base[i]],num_o2[base[i]]);
			mpz_mod(num_r[base[i]], num_r[base[i]],pool[base[i]].modulo);


		}


}

	void RNS_Add_mod_M(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool){

		mpz_t temp1[2*MOD_NUM];
		mpz_t temp2[2*MOD_NUM];

		int i=0;
		for (i=0;i<2*MOD_NUM;i++){
			mpz_init(temp1[i]);
			mpz_init(temp2[i]);
		}


		RNS_add_no(temp1,num_x,num_y,base1,pool);


		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].min_inv_p_mod);
		mpz_mod(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].modulo);

		}



		RNS_add_no(temp1,num_x,num_y,base2,pool);

		Base_extention_B(temp1,temp2,base1,base2,pool);


		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].p_mod);
		mpz_mod(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].modulo);
		}



		RNS_add_no(temp2,temp1,temp2,base2,pool);

		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp1[base2[i]],temp2[base2[i]],pool[base2[i]].Mb_inv_mi[base1[MOD_NUM]]);
		mpz_mod(num_r[base2[i]],temp1[base2[i]],pool[base2[i]].modulo);
		}

		Base_extention_B(temp1, num_r,base2,base1,pool);

		for (i=0;i<2*MOD_NUM;i++){
			mpz_clear(temp1[i]);
			mpz_clear(temp2[i]);
		}

	}
void RNS_Sub_mod_M(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool){

		mpz_t temp1[2*MOD_NUM];
		mpz_t temp2[2*MOD_NUM];

		int i=0;
		for (i=0;i<2*MOD_NUM;i++){
			mpz_init(temp1[i]);
			mpz_init(temp2[i]);
		}


		RNS_sub_no(temp1,num_x,num_y,base1,pool);


		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].min_inv_p_mod);
		mpz_mod(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].modulo);

		}



		RNS_sub_no(temp1,num_x,num_y,base2,pool);



		Base_extention_B(temp1,temp2,base1,base2,pool);

		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].p_mod);
		mpz_mod(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].modulo);
		}



		RNS_add_no(temp2,temp1,temp2,base2,pool);

		for (i=0; i<MOD_NUM;i++){
		mpz_mul(temp1[base2[i]],temp2[base2[i]],pool[base2[i]].Mb_inv_mi[base1[MOD_NUM]]);
		mpz_mod(num_r[base2[i]],temp1[base2[i]],pool[base2[i]].modulo);
		}

		Base_extention_B(temp1, num_r,base2,base1,pool);

		for (i=0;i<2*MOD_NUM;i++){
			mpz_clear(temp1[i]);
			mpz_clear(temp2[i]);
		}

	}
void RNS_Add_mod(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool, mpz_t field_p)
{


	mpz_t num1,num2;
	mpz_init(num1);
	mpz_init(num2);
	rns_to_binary(&num1,num_x,base1,pool);
	rns_to_binary(&num2,num_y,base1,pool);
	mpz_add(num1,num2,num1);
	if (mpz_cmp(field_p,num1)<=0){
		mpz_sub(num1,num1,field_p);
	}
	binary_to_rns(num1, num_r,base1, pool);
	binary_to_rns(num1, num_r,base2, pool);
	mpz_clear(num1);
	mpz_clear(num2);

/*
		mpz_t temp1[2*MOD_NUM];
//		mpz_t rns_p[2*MOD_NUM];
		mpz_t num;
		mpz_init(num);
		int t=0;
		for (t=0;t<2*MOD_NUM;t++){
			mpz_init(temp1[t]);
//			mpz_init(rns_p[i]);
//			mpz_set(rns_p[i],pool[i].p_mod);
		}

		for(t=0;t<MOD_NUM;t++){
			mpz_add(temp1[base1[t]],num_x[base1[t]],num_y[base1[t]]);
			mpz_mod(temp1[base1[t]], temp1[base1[t]],pool[base1[t]].modulo);
			mpz_add(temp1[base2[t]],num_x[base2[t]],num_y[base2[t]]);
			mpz_mod(temp1[base2[t]], temp1[base2[t]],pool[base2[t]].modulo);
			}
		//RNS_add_no(temp1,num_x,num_y,base1,pool);
		//RNS_add_no(temp1,num_x,num_y,base2,pool);
		rns_to_binary(&num,temp1,base1,pool);
		if (mpz_cmp(num, field_p)>=0){
			gmp_printf("%d n>=p %Zd \n    p: %Zd\n",base1[4],num, field_p);
			for(t=0;t<MOD_NUM;t++){
				mpz_sub(temp1[base1[t]],temp1[base1[t]],pool[base1[t]].p_mod);
				mpz_mod(num_r[base1[t]], temp1[base1[t]],pool[base1[t]].modulo);
				mpz_sub(temp1[base2[t]],temp1[base2[t]],pool[base2[t]].p_mod);
				mpz_mod(num_r[base2[t]], temp1[base2[t]],pool[base2[t]].modulo);
			}

			//RNS_sub_no(num_r,temp1,rns_p,base1,pool);
			//RNS_sub_no(num_r,temp1,rns_p,base2,pool);
			rns_to_binary(&num,num_r,base1,pool);
			gmp_printf(" n=   %Zd \n  ++++++++\n",num);
		}

		for (t=0;t<2*MOD_NUM;t++){
			mpz_clear(temp1[t]);
//			mpz_clear(rns_p[i]);
		}
		mpz_clear(num);
*/
}
void RNS_Sub_mod(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool, mpz_t field_p)
{
	mpz_t num1,num2;
	mpz_init(num1);
	mpz_init(num2);
	rns_to_binary(&num1,num_x,base1,pool);
	rns_to_binary(&num2,num_y,base1,pool);
	if (mpz_cmp(num1, num2)<0){
		mpz_add(num1,field_p,num1);
		mpz_sub(num1,num1,num2);
	}else {mpz_sub(num1,num1,num2);}

	binary_to_rns(num1, num_r,base1, pool);
	binary_to_rns(num1, num_r,base2, pool);
	mpz_clear(num1);
	mpz_clear(num2);


/*
		mpz_t temp1[2*MOD_NUM];
		//mpz_t rns_p[2*MOD_NUM];
		mpz_t num;
		mpz_init(num);

		int t=0;
		for (t=0;t<2*MOD_NUM;t++){
			mpz_init(temp1[t]);
//			mpz_init(rns_p[i]);
		}

		for(t=0;t<MOD_NUM;t++){
			mpz_sub(temp1[base1[t]],num_x[base1[t]],num_y[base1[t]]);
			mpz_mod(temp1[base1[t]], temp1[base1[t]],pool[base1[t]].modulo);
			mpz_sub(temp1[base2[t]],num_x[base2[t]],num_y[base2[t]]);
			mpz_mod(temp1[base2[t]], temp1[base2[t]],pool[base2[t]].modulo);
			}

		//RNS_sub_no(temp1,num_x,num_y,base1,pool);
		//RNS_sub_no(temp1,num_x,num_y,base2,pool);
		rns_to_binary(&num,temp1,base1,pool);

		if (mpz_cmp(num, field_p)>=0){
			gmp_printf("%d n>=p %Zd \n   p: %Zd\n",base1[4],num, field_p);

			for(t=0;t<MOD_NUM;t++){
				mpz_add(temp1[base1[t]],temp1[base1[t]],pool[base1[t]].p_mod);
				mpz_mod(num_r[base1[t]], temp1[base1[t]],pool[base1[t]].modulo);
				mpz_add(temp1[base2[t]],temp1[base2[t]],pool[base2[t]].p_mod);
				mpz_mod(num_r[base2[t]], temp1[base2[t]],pool[base2[t]].modulo);
			}

			//RNS_add_no(num_r,temp1,rns_p,base1,pool);
			//RNS_add_no(num_r,temp1,rns_p,base2,pool);

			rns_to_binary(&num,num_r,base1,pool);
		//	mpz_mod(num,num,field_p);
		//	binary_to_rns(num, num_r,base1,pool);
		//	binary_to_rns(num, num_r,base2,pool);
		//	rns_to_binary(&num,num_r,base2,pool);
			gmp_printf(" n=   %Zd \n  ---------\n",num);
		}

		for (t=0;t<2*MOD_NUM;t++){
			mpz_clear(temp1[t]);
//			mpz_clear(rns_p[i]);
		}
		mpz_clear(num);
*/
}

void RNS_Montg_mul(mpz_t num_r[], mpz_t num_x[], mpz_t num_y[], int base1[], int base2[], rns_base_element_data pool){

	mpz_t temp1[2*MOD_NUM];
	mpz_t temp2[2*MOD_NUM];

	int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(temp1[i]);
		mpz_init(temp2[i]);
	}


	RNS_mul(temp1,num_x,num_y,base1,pool);


	for (i=0; i<MOD_NUM;i++){
	mpz_mul(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].min_inv_p_mod);
	mpz_mod(temp1[base1[i]],temp1[base1[i]],pool[base1[i]].modulo);

	}



	RNS_mul(temp1,num_x,num_y,base2,pool);


	Base_extention_B(temp1,temp2,base1,base2,pool);


	for (i=0; i<MOD_NUM;i++){
	mpz_mul(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].p_mod);
	mpz_mod(temp2[base2[i]],temp2[base2[i]],pool[base2[i]].modulo);
	}


	RNS_add_no(temp2,temp1,temp2,base2,pool);


	for (i=0; i<MOD_NUM;i++){
	mpz_mul(temp1[base2[i]],temp2[base2[i]],pool[base2[i]].Mb_inv_mi[base1[MOD_NUM]]);
	mpz_mod(num_r[base2[i]],temp1[base2[i]],pool[base2[i]].modulo);
	}

	Base_extention_B(temp1, num_r,base2,base1,pool);

	for (i=0;i<2*MOD_NUM;i++){
		mpz_clear(temp1[i]);
		mpz_clear(temp2[i]);
	}

}

void EdC_point_add(ec_point_rns *P3, ec_point_rns P1, ec_point_rns P2, int base1[], int base2[], rns_base_element_data pool,mpz_t field_p){


	mpz_t a[2*MOD_NUM];
	mpz_t b[2*MOD_NUM];
	mpz_t c[2*MOD_NUM];
	mpz_t d[2*MOD_NUM];
	mpz_t e[2*MOD_NUM];
	mpz_t f[2*MOD_NUM];
	mpz_t g[2*MOD_NUM];
	mpz_t h[2*MOD_NUM];
	mpz_t q[2*MOD_NUM];
	int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(a[i]);
		mpz_init(b[i]);
		mpz_init(c[i]);
		mpz_init(d[i]);
		mpz_init(e[i]);
		mpz_init(f[i]);
		mpz_init(g[i]);
		mpz_init(h[i]);
		mpz_init(q[i]);
	}


	RNS_Montg_mul(a, P1.z, P2.z,base1,base2,pool);
	RNS_Montg_mul(b, a, a,base1,base2,pool);
	RNS_Montg_mul(c, P1.x, P2.x,base1,base2,pool);
	RNS_Montg_mul(d, P1.y, P2.y, base1,base2,pool);
	RNS_Montg_mul(e, c, d,base1,base2,pool);

	RNS_Add_mod(e,e,e,base1, base2, pool,field_p);		  //to create the 2*C assuming that d=2
	//RNS_Add_mod(e,e,e,base2,pool, field_p);  //to create the 2*C assuming that d=2
	//RNS_Add_mod(e,e,e,base2,pool, field_p);  //to create the 2*C assuming that d=2

	RNS_Sub_mod(f,b,e,base1,base2, pool,field_p);
	//RNS_Sub_mod(f,b,e,base2,pool,field_p);
	//RNS_Sub_mod(f,b,e,base2,pool,field_p);


	RNS_Add_mod(g,b,e,base1,base2,pool,field_p); // b not needed anymore
	//RNS_Add_mod(g,b,e,base2,pool,field_p); // Tests, Just ignore
	//RNS_Add_mod(g,b,e,base2,pool,field_p); // Tests, Just ignore

	RNS_Add_mod(b,c,d,base1,base2,pool,field_p);
	//RNS_Add_mod(b,c,d,base2,pool,field_p); //Tests, Just ignore
	//RNS_Add_mod(b,c,d,base2,pool,field_p); //Tests, Just ignore

	RNS_Sub_mod(c,d,c,base1,base2,pool,field_p);
	//RNS_Sub_mod(c,d,c,base2,pool,field_p); //Tests, Just ignore
	//RNS_Sub_mod(c,d,c,base2,pool,field_p); //Tests, Just ignore

	RNS_Add_mod(d,P1.x,P1.y,base1,base2,pool,field_p);
	//RNS_Add_mod(d,P1.x,P1.y,base2,pool,field_p); //Tests, Just ignore
	//RNS_Add_mod(d,P1.x,P1.y,base2,pool,field_p); //Tests, Just ignore

	RNS_Add_mod(h,P2.x,P2.y,base1,base2,pool,field_p);
	//RNS_Add_mod(h,P2.x,P2.y,base2,pool,field_p); //Tests, Just ignore
	//RNS_Add_mod(h,P2.x,P2.y,base2,pool,field_p); //Tests, Just ignore

	RNS_Montg_mul(q, h, d, base1,base2,pool);

	RNS_Sub_mod(h,q,b,base1,base2,pool,field_p);
	//RNS_Sub_mod(h,q,b,base2,pool,field_p); //Tests, Just ignore
	//RNS_Sub_mod(h,q,b,base2,pool,field_p); //Tests, Just ignore

	RNS_Montg_mul(q, h, f,base1,base2,pool);

	RNS_Montg_mul(P3->x,q, a,base1,base2,pool);

	RNS_Montg_mul(q,g, c, base1,base2,pool);
	RNS_Montg_mul(P3->y,q, a,base1,base2,pool);

	RNS_Montg_mul(P3->z,f, g,base1,base2,pool); //assuming c =1

	for (i=0;i<2*MOD_NUM;i++){
		mpz_clear(a[i]);
		mpz_clear(b[i]);
		mpz_clear(c[i]);
		mpz_clear(d[i]);
		mpz_clear(e[i]);
		mpz_clear(f[i]);
		mpz_clear(g[i]);
		mpz_clear(h[i]);
		mpz_clear(q[i]);
	}


}

void EdC_point_double(ec_point_rns *P3, ec_point_rns P1, int base1[], int base2[], rns_base_element_data pool, mpz_t field_p){

	mpz_t a[2*MOD_NUM];
	mpz_t b[2*MOD_NUM];
	mpz_t c[2*MOD_NUM];
	mpz_t d[2*MOD_NUM];
	mpz_t e[2*MOD_NUM];
	mpz_t f[2*MOD_NUM];
	mpz_t g[2*MOD_NUM];

	int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(a[i]);
		mpz_init(b[i]);
		mpz_init(c[i]);
		mpz_init(d[i]);
		mpz_init(e[i]);
		mpz_init(f[i]);
		mpz_init(g[i]);

	}


	RNS_Add_mod(a,P1.x,P1.y,base1,base2,pool,field_p);
	//RNS_Add_mod_M(a,P1.x,P1.y,base1,base2,pool); //Tests, Just ignore
	//RNS_add(a,P1.x,P1.y,base1,pool); //Tests, Just ignore
	//RNS_add(a,P1.x,P1.y,base2,pool); //Tests, Just ignore

	RNS_Montg_mul(b,a, a,base1,base2,pool);

	RNS_Montg_mul(c,P1.x,P1.x,base1,base2,pool);
	RNS_Montg_mul(d,P1.y, P1.y, base1,base2,pool);

	RNS_Montg_mul(a,P1.z,P1.z,base1,base2,pool);

	RNS_Add_mod(e,c,d,base1,base2,pool,field_p);
	//RNS_Add_mod_M(e,c,d,base1,base2,pool); //Tests, Just ignore
	//RNS_add(e,c,d,base1,pool); //Tests, Just ignore
	//RNS_add(e,c,d,base2,pool); //Tests, Just ignore

	RNS_Add_mod(a,a,a,base1,base2,pool,field_p);
	//RNS_Add_mod_M(a,a,a,base1,base2,pool); //Tests, Just ignore
	//RNS_add(a,a,a,base1,pool);  //Tests, Just ignore
	//RNS_add(a,a,a,base2,pool);  //Tests, Just ignore

	RNS_Sub_mod(f,e,a,base1,base2,pool,field_p);
	//RNS_Sub_mod_M(f,e,a,base1,base2,pool); //Tests, Just ignore
	//RNS_sub(f,e,a,base1,pool); 		//Tests, Just ignore
	//RNS_sub(f,e,a,base2,pool); 		//Tests, Just ignore

	RNS_Sub_mod(b,b,e,base1,base2,pool,field_p);
	//RNS_Sub_mod_M(b,b,e,base1,base2,pool); //Tests, Just ignore
	//RNS_sub(b,b,e,base1,pool);   //Tests, Just ignore
	//RNS_sub(b,b,e,base2,pool);   //Tests, Just ignore

	RNS_Montg_mul(P3->x,b,f,base1,base2,pool);


	RNS_Sub_mod(g,c,d,base1,base2,pool,field_p);
	//RNS_Sub_mod_M(g,c,d,base1,base2,pool);  //Tests, Just ignore
	//RNS_sub(g,c,d,base1,pool); //Tests, Just ignore
	//RNS_sub(g,c,d,base2,pool); //Tests, Just ignore



	RNS_Montg_mul(P3->y,g,e,base1,base2,pool);

	RNS_Montg_mul(P3->z, e, f, base1,base2,pool);

	for (i=0;i<2*MOD_NUM;i++){
		mpz_clear(a[i]);
		mpz_clear(b[i]);
		mpz_clear(c[i]);
		mpz_clear(d[i]);
		mpz_clear(e[i]);
		mpz_clear(f[i]);

	}

}


void Rand_EC_point(ec_point *Q, ec_point base_pointP, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p){
	struct timeval  tv1;
	int n=192;
	mpz_t scalar;
	  gmp_randstate_t exp_seed;
	gmp_randinit_default(exp_seed);
	 gettimeofday(&tv1,NULL);
	 gmp_randseed_ui(exp_seed,tv1.tv_usec);
	 mpz_urandomb(scalar, exp_seed, (mp_bitcnt_t ) n);
	 mpz_out_str(stdout, 10, scalar);

	 EC_scalar_mul_no_Base_Perm_no_R_point(Q,base_pointP,scalar,pool,base_rand_index,modp,field_p);

}

void EC_scalar_mul(ec_point *Q, ec_point P, ec_point R, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p){

	int Bn[MOD_NUM+1];
	int Bnn[MOD_NUM+1];

	int init_Bn[MOD_NUM+1];
	int init_Bnn[MOD_NUM+1];

	int Bn_new[MOD_NUM+1];
	int Bnn_new[MOD_NUM+1];

	ec_point_rns R0=rns_EC_point_init();
	ec_point_rns R1=rns_EC_point_init();
	ec_point_rns R2=rns_EC_point_init();



	mpz_t mmodp_B[2*MOD_NUM];
	mpz_t val_rns[2*MOD_NUM];
	mpz_t val_rns_1[2*MOD_NUM];


	int i=0;
	for (i=0; i<2*MOD_NUM;i++){


		mpz_init(mmodp_B[i]);
		mpz_init(val_rns[i]);
		mpz_init(val_rns_1[i]);
	}

	struct timeval  tv1;
		int n=192;
	mpz_t permut_index;
	mpz_init(permut_index);

	gmp_randstate_t state;
	gettimeofday(&tv1,NULL);
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,tv1.tv_usec);

	mpz_set_str(val_rns[0],"70",10);
	mpz_urandomm(permut_index, state, val_rns[0]);

	int pem_in=(int)mpz_get_ui(permut_index);
	int init_permut_index=pem_in;
	int rand_checker[70];
	for (i=0;i<70;i++){
		rand_checker[i]=0;
	}
	//int init_permut_index=0;
	//int pem_in=0;
	for (i=0;i<MOD_NUM;i++){
		init_Bn[i]=base_rand_index[init_permut_index].moduli[i];
		init_Bnn[i]=base_rand_index[init_permut_index].m_nn_nums[i];
		Bn[i]=base_rand_index[init_permut_index].moduli[i];
		Bnn[i]=base_rand_index[init_permut_index].m_nn_nums[i];
		}
		init_Bn[MOD_NUM]=base_rand_index[init_permut_index].base_index;
		init_Bnn[MOD_NUM]=69-base_rand_index[init_permut_index].base_index;
		Bn[MOD_NUM]=base_rand_index[init_permut_index].base_index;
		Bnn[MOD_NUM]=69-base_rand_index[init_permut_index].base_index;

		for (i=0; i<2*MOD_NUM;i++){
			mpz_set_str(val_rns[i],"0",10);
			mpz_set_str(val_rns_1[i],"1",10);
		}
		for (i=0;i<MOD_NUM;i++){
			mpz_set(mmodp_B[Bn[i]],modp.m[Bn[i]]);
			mpz_set(mmodp_B[Bnn[i]],modp.m[Bnn[i]]);
		}


		//base point in rns format

		ec_point V=EC_point_init();

		binary_to_rns(P.x, R1.x, Bn, pool);
		binary_to_rns(P.x, R1.x, Bnn, pool);
		binary_to_rns(P.y, R1.y, Bn, pool);
		binary_to_rns(P.y, R1.y, Bnn, pool);
		binary_to_rns(P.z, R1.z, Bn, pool);
		binary_to_rns(P.z, R1.z, Bnn, pool);


		binary_to_rns(R.x, R0.x, Bn, pool);
		binary_to_rns(R.x, R0.x, Bnn, pool);
		binary_to_rns(R.y, R0.y, Bn, pool);
		binary_to_rns(R.y, R0.y, Bnn, pool);
		binary_to_rns(R.z, R0.z, Bn, pool);
		binary_to_rns(R.z, R0.z, Bnn, pool);

		mpz_sub(R.x,val_rns[0],R.x);
		//mpz_mod(R.x,R.x,field_p);

		binary_to_rns(R.x, R2.x, Bn, pool);
		binary_to_rns(R.x, R2.x, Bnn, pool);




//conversion to Montgomery format


		RNS_Montg_mul(R1.x, R1.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.y, R1.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.z, R1.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R0.x, R0.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.y, R0.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.z, R0.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R2.x, R2.x, mmodp_B,Bnn,Bn,pool);



		for (i=0; i<2*MOD_NUM;i++){

		mpz_set(R2.y[i],R0.y[i]);
		mpz_set(R2.z[i],R0.z[i]);
		}



		EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
		gmp_printf("----\n");
//----------------------------------
#ifdef DEB_MONT
		ec_point_rns output_point_rns=rns_EC_point_init();
		EdC_point_add(&output_point_rns, R2, R0, Bn, Bnn,pool,field_p);
		rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
		rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
		rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
		mpz_mod(V.x,V.x,field_p);
		mpz_mod(V.y,V.y,field_p);
		mpz_mod(V.z,V.z,field_p);
		gmp_printf(" check Q= %Zd, %Zd, %Zd\n",V.x,V.y,V.z);
//----------------------------------
#endif

	for (i=EC_Field_P_bit_lenght-1;i>=0;i--){
		int j;
		EdC_point_double(&R2, R2, init_Bn, init_Bnn,pool,field_p);  // TODO problem

		//-----------Choose random permutation --------



			mpz_set_str(val_rns[0],"70",10);
			mpz_urandomm(permut_index, state, val_rns[0]);
		    pem_in=(int)mpz_get_ui(permut_index);

			mpz_set_str(val_rns[0],"0",10);
			rand_checker[pem_in]++;

			for (j=0; j<MOD_NUM;j++){
				Bn_new[j]=base_rand_index[pem_in].moduli[j];
				Bnn_new[j]=base_rand_index[pem_in].m_nn_nums[j];
			}
			Bn_new[MOD_NUM]=base_rand_index[pem_in].base_index;
			Bnn_new[MOD_NUM]=69-base_rand_index[pem_in].base_index;



			RNS_Montg_mul(R1.x,R1.x, mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.x, R1.x,val_rns_1, Bn,Bnn,pool);

			RNS_Montg_mul(R1.y,R1.y,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.y, R1.y, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R1.z,R1.z,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.z, R1.z, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R0.x, R0.x,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.x, R0.x, val_rns_1,Bn,Bnn,pool);


			RNS_Montg_mul(R0.y,R0.y,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.y, R0.y, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R0.z,R0.z,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.z, R0.z,val_rns_1,Bn,Bnn,pool);



			for (j=0;j<=MOD_NUM;j++){
				Bn[j]=Bn_new[j];
				Bnn[j]=Bnn_new[j];
			}



#ifdef DEB_MONT
			RNS_Montg_mul(output_point_rns.x, R0.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R0.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R0.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("T0= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);

			RNS_Montg_mul(output_point_rns.x, R1.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R1.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R1.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("T1= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);
			printf("K = point_add(P[0], P[1], P[2], T0[0], T0[1], T0[2], d, Fp)\n");
			printf("ss=10011%d\n",i);
			printf("sss=%d\n",Bn[4]);
			printf("print(gmpy2.f_mod(K[0] * gmpy2.invert(K[2], Fp), Fp), gmpy2.f_mod(K[1] * gmpy2.invert(K[2], Fp), Fp))\n");
			printf("print(gmpy2.f_mod(T1[0] * gmpy2.invert(T1[2], Fp), Fp), gmpy2.f_mod(T1[1] * gmpy2.invert(T1[2], Fp), Fp))\n");
			printf("print(ss,sss)\n");
#endif



		if (mpz_tstbit(k,i)==1){

			EdC_point_add(&R0, R1, R0, Bn, Bnn,pool,field_p);
#ifdef DEB_MONT
			RNS_Montg_mul(output_point_rns.x, R1.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R1.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R1.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("R1= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);
			printf("print(R1[0],R1[1],R1[2])\n");
#endif

			EdC_point_double(&R1, R1, Bn, Bnn,pool,field_p);
#ifdef DEB_MONT
			RNS_Montg_mul(output_point_rns.x, R0.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R0.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R0.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("R0= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);

			RNS_Montg_mul(output_point_rns.x, R1.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R1.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R1.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("R1= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);
			printf("K = point_add(P[0], P[1], P[2], R0[0], R0[1], R0[2], d, Fp)\n");
			printf("ss=111%d\n",i);
			printf("sss=%d\n",Bn[4]);
			printf("print(gmpy2.f_mod(K[0] * gmpy2.invert(K[2], Fp), Fp), gmpy2.f_mod(K[1] * gmpy2.invert(K[2], Fp), Fp))\n");
			printf("print(gmpy2.f_mod(R1[0] * gmpy2.invert(R1[2], Fp), Fp), gmpy2.f_mod(R1[1] * gmpy2.invert(R1[2], Fp), Fp))\n");
			printf("print(ss,sss)\n");
#endif
		}
		else {

			EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
			EdC_point_double(&R0, R0, Bn, Bnn,pool,field_p);
#ifdef DEB_MONT
			RNS_Montg_mul(output_point_rns.x, R0.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R0.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R0.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("R0= %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);

			RNS_Montg_mul(output_point_rns.x, R1.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.y, R1.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
			RNS_Montg_mul(output_point_rns.z, R1.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


			rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
			rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
			rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
			gmp_printf("R1 = %Zd, %Zd, %Zd  \n",V.x, V.y, V.z);
			printf("K = point_add(P[0], P[1], P[2], R0[0], R0[1], R0[2], d, Fp)\n");
			printf("ss=101%d\n",i);
			printf("sss=%d\n",Bn[4]);
			printf("print(gmpy2.f_mod(K[0] * gmpy2.invert(K[2], Fp), Fp), gmpy2.f_mod(K[1] * gmpy2.invert(K[2], Fp), Fp))\n");
			printf("print(gmpy2.f_mod(R1[0] * gmpy2.invert(R1[2], Fp), Fp), gmpy2.f_mod(R1[1] * gmpy2.invert(R1[2], Fp), Fp))\n");
			printf("print(ss,sss)\n");
#endif
		}
	}


	ec_point_rns tempP=rns_EC_point_init();

	 binary_to_rns(P.x, tempP.x, Bn, pool);
	 binary_to_rns(P.x, tempP.x, Bnn, pool);
	 binary_to_rns(P.y, tempP.y, Bn, pool);
	 binary_to_rns(P.y, tempP.y, Bnn, pool);
	 binary_to_rns(P.z, tempP.z, Bn, pool);
	 binary_to_rns(P.z, tempP.z, Bnn, pool);
#ifdef DEB_MONT
		rns_to_binary(&V.x, tempP.x,Bn,pool);
		rns_to_binary(&V.y, tempP.y,Bn,pool);
		rns_to_binary(&V.z, tempP.z,Bn,pool);
		gmp_printf("P[%d] Q=%Zd, %Zd, %Zd  \n",Bn[4],V.x, V.y, V.z);
#endif
	RNS_Montg_mul(tempP.x, tempP.x, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.y, tempP.y, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.z, tempP.z, mmodp_B,Bnn,Bn,pool);

	 EdC_point_add(&tempP, tempP, R0, Bn, Bnn,pool,field_p);
#ifdef DEB_MONT
		RNS_Montg_mul(output_point_rns.x, R0.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
 		RNS_Montg_mul(output_point_rns.y, R0.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
  		RNS_Montg_mul(output_point_rns.z, R0.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only


		rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
		rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
		rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
		gmp_printf("R0[%d] Q=%Zd, %Zd, %Zd  \n",Bn[4],V.x, V.y, V.z);

		RNS_Montg_mul(output_point_rns.x, tempP.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
  		RNS_Montg_mul(output_point_rns.y, tempP.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
 		RNS_Montg_mul(output_point_rns.z, tempP.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only

  		rns_to_binary(&V.x, output_point_rns.x,Bn,B_pool);
		rns_to_binary(&V.y, output_point_rns.y,Bn,B_pool);
		rns_to_binary(&V.z, output_point_rns.z,Bn,B_pool);
		gmp_printf("R0+P Q=%Zd, %Zd, %Zd  \n",V.x, V.y, V.z);


		RNS_Montg_mul(output_point_rns.x, R1.x, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
  		RNS_Montg_mul(output_point_rns.y, R1.y, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only
 		RNS_Montg_mul(output_point_rns.z, R1.z, val_rns_1,Bn,Bnn, B_pool);	 // not needed for testing only

		rns_to_binary(&V.x, output_point_rns.x,Bn,pool);
		rns_to_binary(&V.y, output_point_rns.y,Bn,pool);
		rns_to_binary(&V.z, output_point_rns.z,Bn,pool);
		gmp_printf("R1 Q=%Zd, %Zd, %Zd  \n",V.x, V.y, V.z);

#endif
	rns_to_binary(&V.x, R1.x,Bn,pool);
	mpz_sub(V.x,val_rns[0],V.x);
	mpz_mod(V.x,V.x,field_p);
	binary_to_rns(V.x,R1.x,Bn,B_pool);
	binary_to_rns(V.x,R1.x,Bnn,B_pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bn,pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bnn,pool);

	EdC_point_add(&tempP, tempP, R1, Bn, Bnn,pool,field_p);

#ifdef DEB_MONT
	rns_to_binary(&V.x, tempP.x,Bn,pool);
	rns_to_binary(&V.y, tempP.y,Bn,pool);
	rns_to_binary(&V.z, tempP.z,Bn,pool);
	gmp_printf("V= %Zd,  %Zd,  %Zd \n",V.x, V.y, V.z);
	//gmp_printf("C (%Zd)  %Zd, (%Zd) %Zd, (%Zd) %Zd \n", tempP.x, tempP.xn, tempP.y, tempP.yn, tempP.z, tempP.zn);

	for (i=0;i<70;i++){
		printf("%d ",rand_checker[i]);
	}
#endif
	if((mpz_cmp(tempP.x[0],val_rns[0])==0) &&
		(mpz_cmp(tempP.x[1],val_rns[1])==0) &&
		(mpz_cmp(tempP.x[2],val_rns[2])==0) &&
		(mpz_cmp(tempP.x[3],val_rns[3])==0) &&
		(mpz_cmp(tempP.x[4],val_rns[4])==0) &&
		(mpz_cmp(tempP.x[5],val_rns[5])==0) &&
		(mpz_cmp(tempP.x[6],val_rns[6])==0) &&
		(mpz_cmp(tempP.x[7],val_rns[7])==0))  //&&

		{

		for (i=0; i<MOD_NUM;i++){
			Bn_new[i]=base_index[init_permut_index].moduli[i];
			Bnn_new[i]=base_index[init_permut_index].m_nn_nums[i];
		}
		Bn_new[MOD_NUM]=base_index[init_permut_index].base_index;
		Bnn_new[MOD_NUM]=69-base_index[init_permut_index].base_index;


		    RNS_Montg_mul(R2.x, R2.x,mmodp_B,Bnn,Bn,pool);
		    RNS_Montg_mul(R2.x, R2.x,val_rns_1,Bn_new,Bnn_new,pool);

		    RNS_Montg_mul(R2.y,R2.y,mmodp_B,Bnn,Bn,pool);
			RNS_Montg_mul(R2.y, R2.y,val_rns_1,Bn_new,Bnn_new,pool);

			RNS_Montg_mul(R2.z,R2.z,mmodp_B,Bnn,Bn,pool);
			RNS_Montg_mul(R2.z, R2.z,val_rns_1,Bn_new,Bnn_new,pool);

			EdC_point_add(&tempP, R2, R0, Bn, Bnn, pool,field_p);

			RNS_Montg_mul(tempP.x, tempP.x,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.y, tempP.y,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.z, tempP.z,val_rns_1,Bn,Bnn,pool);

	}
	else
	{ for (i=0; i<2*MOD_NUM;i++){
		mpz_set_str(tempP.x[i],"0",10);
		mpz_set_str(tempP.y[i],"1",10);
		mpz_set_str(tempP.z[i],"1",10);
	 	 }

	printf("faulty\n");

	}
rns_to_binary( &Q->x,tempP.x,Bn,pool);
rns_to_binary(&Q->y,tempP.y,Bn,pool);
rns_to_binary(&Q->z,tempP.z,Bn,pool);

ec_point_rns_clear(tempP);
ec_point_rns_clear(R0);
ec_point_rns_clear(R1);
ec_point_rns_clear(R2);
ec_point_clear(V);

for (i=0;i<2*MOD_NUM;i++){
	mpz_clear(val_rns[i]);
	mpz_clear(val_rns_1[i]);
	mpz_clear(mmodp_B[i]);

}
mpz_clear(permut_index);
	gmp_randclear(state);
}
//#endif

void rns_num_clear(rns_num t){
	mpz_clear(t.m0);
	mpz_clear(t.m1);
	mpz_clear(t.m2);
	mpz_clear(t.m3);
}

void ec_point_clear(ec_point t){
	mpz_clear(t.x);
	mpz_clear(t.y);
	mpz_clear(t.z);
}

void ec_point_rns_clear(ec_point_rns t){
	int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_clear(t.x[i]);
		mpz_clear(t.y[i]);
		mpz_clear(t.z[i]);
	}

}

void EC_scalar_mul_no_Base_Perm(ec_point *Q, ec_point P, ec_point R, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p){

	int Bn[MOD_NUM+1];
	int Bnn[MOD_NUM+1];


	ec_point_rns R0=rns_EC_point_init();
	ec_point_rns R1=rns_EC_point_init();
	ec_point_rns R2=rns_EC_point_init();



	mpz_t mmodp_B[2*MOD_NUM];
	mpz_t val_rns[2*MOD_NUM];
	mpz_t val_rns_1[2*MOD_NUM];


	int i=0;
	for (i=0; i<2*MOD_NUM;i++){


		mpz_init(mmodp_B[i]);
		mpz_init(val_rns[i]);
		mpz_init(val_rns_1[i]);
	}

	int  permut_index=1;  //choise of RNS bases to be used in all scalar multiplication
										// Currently RNS bases choice 1


	for (i=0;i<MOD_NUM;i++){
		Bn[i]=base_rand_index[permut_index].moduli[i];
		Bnn[i]=base_rand_index[permut_index].m_nn_nums[i];
		}
		Bn[MOD_NUM]=base_rand_index[permut_index].base_index;
		Bnn[MOD_NUM]=69-base_rand_index[permut_index].base_index;

		for (i=0; i<2*MOD_NUM;i++){
			mpz_set_str(val_rns[i],"0",10);
			mpz_set_str(val_rns_1[i],"1",10);
		}
		for (i=0;i<MOD_NUM;i++){
			mpz_set(mmodp_B[Bn[i]],modp.m[Bn[i]]);
			mpz_set(mmodp_B[Bnn[i]],modp.m[Bnn[i]]);
		}


		//base point in rns format

		ec_point V=EC_point_init();

		ec_point_rns tempP=rns_EC_point_init();

			 binary_to_rns(P.x, tempP.x, Bn, pool);
			 binary_to_rns(P.x, tempP.x, Bnn, pool);
			 binary_to_rns(P.y, tempP.y, Bn, pool);
			 binary_to_rns(P.y, tempP.y, Bnn, pool);
			 binary_to_rns(P.z, tempP.z, Bn, pool);
			 binary_to_rns(P.z, tempP.z, Bnn, pool);

			 for (i=0;i<2*MOD_NUM;i++){
			 mpz_set(R1.x[i],tempP.x[i]);
			 mpz_set(R1.y[i],tempP.y[i]);
			 mpz_set(R1.z[i],tempP.z[i]);
		}



		binary_to_rns(R.x, R0.x, Bn, pool);
		binary_to_rns(R.x, R0.x, Bnn, pool);
		binary_to_rns(R.y, R0.y, Bn, pool);
		binary_to_rns(R.y, R0.y, Bnn, pool);
		binary_to_rns(R.z, R0.z, Bn, pool);
		binary_to_rns(R.z, R0.z, Bnn, pool);

		mpz_sub(R.x,val_rns[0],R.x);
		//mpz_mod(R.x,R.x,field_p);

		binary_to_rns(R.x, R2.x, Bn, pool);
		binary_to_rns(R.x, R2.x, Bnn, pool);




//conversion to Montgomery format


		RNS_Montg_mul(R1.x, R1.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.y, R1.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.z, R1.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R0.x, R0.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.y, R0.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.z, R0.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R2.x, R2.x, mmodp_B,Bnn,Bn,pool);



		for (i=0; i<2*MOD_NUM;i++){

		mpz_set(R2.y[i],R0.y[i]);
		mpz_set(R2.z[i],R0.z[i]);
		}



		EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
		gmp_printf("----\n");


	for (i=EC_Field_P_bit_lenght-1;i>=0;i--){
		int j;
		EdC_point_double(&R2, R2, Bn, Bnn,pool,field_p);


		if (mpz_tstbit(k,i)==1){

			EdC_point_add(&R0, R1, R0, Bn, Bnn,pool,field_p);

			EdC_point_double(&R1, R1, Bn, Bnn,pool,field_p);

		}
		else {

			EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
			EdC_point_double(&R0, R0, Bn, Bnn,pool,field_p);
		}
	}




	RNS_Montg_mul(tempP.x, tempP.x, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.y, tempP.y, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.z, tempP.z, mmodp_B,Bnn,Bn,pool);

	 EdC_point_add(&tempP, tempP, R0, Bn, Bnn,pool,field_p);

	rns_to_binary(&V.x, R1.x,Bn,pool);
	mpz_sub(V.x,val_rns[0],V.x);
	mpz_mod(V.x,V.x,field_p);
	binary_to_rns(V.x,R1.x,Bn,B_pool);
	binary_to_rns(V.x,R1.x,Bnn,B_pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bn,pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bnn,pool);

	EdC_point_add(&tempP, tempP, R1, Bn, Bnn,pool,field_p);


	if((mpz_cmp(tempP.x[0],val_rns[0])==0) &&
		(mpz_cmp(tempP.x[1],val_rns[1])==0) &&
		(mpz_cmp(tempP.x[2],val_rns[2])==0) &&
		(mpz_cmp(tempP.x[3],val_rns[3])==0) &&
		(mpz_cmp(tempP.x[4],val_rns[4])==0) &&
		(mpz_cmp(tempP.x[5],val_rns[5])==0) &&
		(mpz_cmp(tempP.x[6],val_rns[6])==0) &&
		(mpz_cmp(tempP.x[7],val_rns[7])==0))  //&&

		{


			EdC_point_add(&tempP, R2, R0, Bn, Bnn, pool,field_p);

			RNS_Montg_mul(tempP.x, tempP.x,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.y, tempP.y,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.z, tempP.z,val_rns_1,Bn,Bnn,pool);

	}
	else
	{ for (i=0; i<2*MOD_NUM;i++){
		mpz_set_str(tempP.x[i],"0",10);
		mpz_set_str(tempP.y[i],"1",10);
		mpz_set_str(tempP.z[i],"1",10);
	 	 }

	printf("faulty\n");

	}
rns_to_binary( &Q->x,tempP.x,Bn,pool);
rns_to_binary(&Q->y,tempP.y,Bn,pool);
rns_to_binary(&Q->z,tempP.z,Bn,pool);

ec_point_rns_clear(tempP);
ec_point_rns_clear(R0);
ec_point_rns_clear(R1);
ec_point_rns_clear(R2);
ec_point_clear(V);

for (i=0;i<2*MOD_NUM;i++){
	mpz_clear(val_rns[i]);
	mpz_clear(val_rns_1[i]);
	mpz_clear(mmodp_B[i]);

}

}


void EC_scalar_mul_no_R_point(ec_point *Q, ec_point P, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p){

	int Bn[MOD_NUM+1];
	int Bnn[MOD_NUM+1];

	int init_Bn[MOD_NUM+1];
	int init_Bnn[MOD_NUM+1];

	int Bn_new[MOD_NUM+1];
	int Bnn_new[MOD_NUM+1];

	ec_point_rns R0=rns_EC_point_init();
	ec_point_rns R1=rns_EC_point_init();
	//ec_point_rns R2=rns_EC_point_init();



	mpz_t mmodp_B[2*MOD_NUM];
	mpz_t val_rns[2*MOD_NUM];
	mpz_t val_rns_1[2*MOD_NUM];


	int i=0;
	for (i=0; i<2*MOD_NUM;i++){


		mpz_init(mmodp_B[i]);
		mpz_init(val_rns[i]);
		mpz_init(val_rns_1[i]);
	}

/*	mpz_t permut_index;
	mpz_init(permut_index);

	gmp_randstate_t state;
	gmp_randinit_mt(state);


	mpz_set_str(val_rns[0],"70",10);
	mpz_urandomm(permut_index, state, val_rns[0]);

	int pem_in=(int)mpz_get_ui(permut_index);
	int init_permut_index=pem_in;
	int rand_checker[70];
	for (i=0;i<70;i++){
		rand_checker[i]=0;
	}*/


	struct timeval  tv1;
		int n=192;
	mpz_t permut_index;
	mpz_init(permut_index);

	gmp_randstate_t state;
	gettimeofday(&tv1,NULL);
	gmp_randinit_mt(state);
	gmp_randseed_ui(state,tv1.tv_usec);

	mpz_set_str(val_rns[0],"70",10);
	mpz_urandomm(permut_index, state, val_rns[0]);

	int pem_in=(int)mpz_get_ui(permut_index);
	int init_permut_index=pem_in;
	int rand_checker[70];
	for (i=0;i<70;i++){
		rand_checker[i]=0;
	}
	//int init_permut_index=0;
	//int pem_in=0;
	for (i=0;i<MOD_NUM;i++){
		init_Bn[i]=base_rand_index[init_permut_index].moduli[i];
		init_Bnn[i]=base_rand_index[init_permut_index].m_nn_nums[i];
		Bn[i]=base_rand_index[init_permut_index].moduli[i];
		Bnn[i]=base_rand_index[init_permut_index].m_nn_nums[i];
		}
		init_Bn[MOD_NUM]=base_rand_index[init_permut_index].base_index;
		init_Bnn[MOD_NUM]=69-base_rand_index[init_permut_index].base_index;
		Bn[MOD_NUM]=base_rand_index[init_permut_index].base_index;
		Bnn[MOD_NUM]=69-base_rand_index[init_permut_index].base_index;

		for (i=0; i<2*MOD_NUM;i++){
			mpz_set_str(val_rns[i],"0",10);
			mpz_set_str(val_rns_1[i],"1",10);
		}
		for (i=0;i<MOD_NUM;i++){
			mpz_set(mmodp_B[Bn[i]],modp.m[Bn[i]]);
			mpz_set(mmodp_B[Bnn[i]],modp.m[Bnn[i]]);
		}


		//base point in rns format

		ec_point V=EC_point_init();

		binary_to_rns(P.x, R1.x, Bn, pool);
		binary_to_rns(P.x, R1.x, Bnn, pool);
		binary_to_rns(P.y, R1.y, Bn, pool);
		binary_to_rns(P.y, R1.y, Bnn, pool);
		binary_to_rns(P.z, R1.z, Bn, pool);
		binary_to_rns(P.z, R1.z, Bnn, pool);


		for (i=0; i<2*MOD_NUM;i++){
		mpz_set_str(R0.x[i],"0",10);
		mpz_set_str(R0.y[i],"1",10);
		mpz_set_str(R0.z[i],"1",10);
		}


//conversion to Montgomery format


		RNS_Montg_mul(R1.x, R1.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.y, R1.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.z, R1.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R0.x, R0.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.y, R0.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.z, R0.z, mmodp_B,Bnn,Bn,pool);


		gmp_printf("----\n");

	for (i=EC_Field_P_bit_lenght-1;i>=0;i--){
		int j;


		//-----------Choose random permutation --------



			mpz_set_str(val_rns[0],"70",10);
			mpz_urandomm(permut_index, state, val_rns[0]);
		    pem_in=(int)mpz_get_ui(permut_index);

			mpz_set_str(val_rns[0],"0",10);
			rand_checker[pem_in]++;

			for (j=0; j<MOD_NUM;j++){
				Bn_new[j]=base_rand_index[pem_in].moduli[j];
				Bnn_new[j]=base_rand_index[pem_in].m_nn_nums[j];
			}
			Bn_new[MOD_NUM]=base_rand_index[pem_in].base_index;
			Bnn_new[MOD_NUM]=69-base_rand_index[pem_in].base_index;



			RNS_Montg_mul(R1.x,R1.x, mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.x, R1.x,val_rns_1, Bn,Bnn,pool);

			RNS_Montg_mul(R1.y,R1.y,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.y, R1.y, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R1.z,R1.z,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R1.z, R1.z, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R0.x, R0.x,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.x, R0.x, val_rns_1,Bn,Bnn,pool);


			RNS_Montg_mul(R0.y,R0.y,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.y, R0.y, val_rns_1,Bn,Bnn,pool);

			RNS_Montg_mul(R0.z,R0.z,mmodp_B,Bnn_new,Bn_new,pool);
			RNS_Montg_mul(R0.z, R0.z,val_rns_1,Bn,Bnn,pool);



			for (j=0;j<=MOD_NUM;j++){
				Bn[j]=Bn_new[j];
				Bnn[j]=Bnn_new[j];
			}




		if (mpz_tstbit(k,i)==1){

			EdC_point_add(&R0, R1, R0, Bn, Bnn,pool,field_p);

			EdC_point_double(&R1, R1, Bn, Bnn,pool,field_p);

		}
		else {

			EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
			EdC_point_double(&R0, R0, Bn, Bnn,pool,field_p);

		}
	}


	ec_point_rns tempP=rns_EC_point_init();

	 binary_to_rns(P.x, tempP.x, Bn, pool);
	 binary_to_rns(P.x, tempP.x, Bnn, pool);
	 binary_to_rns(P.y, tempP.y, Bn, pool);
	 binary_to_rns(P.y, tempP.y, Bnn, pool);
	 binary_to_rns(P.z, tempP.z, Bn, pool);
	 binary_to_rns(P.z, tempP.z, Bnn, pool);

	RNS_Montg_mul(tempP.x, tempP.x, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.y, tempP.y, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.z, tempP.z, mmodp_B,Bnn,Bn,pool);

	 EdC_point_add(&tempP, tempP, R0, Bn, Bnn,pool,field_p);

	rns_to_binary(&V.x, R1.x,Bn,pool);
	mpz_sub(V.x,val_rns[0],V.x);
	mpz_mod(V.x,V.x,field_p);
	binary_to_rns(V.x,R1.x,Bn,B_pool);
	binary_to_rns(V.x,R1.x,Bnn,B_pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bn,pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bnn,pool);

	EdC_point_add(&tempP, tempP, R1, Bn, Bnn,pool,field_p);


	if((mpz_cmp(tempP.x[0],val_rns[0])==0) &&
		(mpz_cmp(tempP.x[1],val_rns[1])==0) &&
		(mpz_cmp(tempP.x[2],val_rns[2])==0) &&
		(mpz_cmp(tempP.x[3],val_rns[3])==0) &&
		(mpz_cmp(tempP.x[4],val_rns[4])==0) &&
		(mpz_cmp(tempP.x[5],val_rns[5])==0) &&
		(mpz_cmp(tempP.x[6],val_rns[6])==0) &&
		(mpz_cmp(tempP.x[7],val_rns[7])==0))  //&&

		{

			RNS_Montg_mul(tempP.x, R0.x,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.y, R0.y,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.z, R0.z,val_rns_1,Bn,Bnn,pool);

	}
	else
	{ for (i=0; i<2*MOD_NUM;i++){
		mpz_set_str(tempP.x[i],"0",10);
		mpz_set_str(tempP.y[i],"1",10);
		mpz_set_str(tempP.z[i],"1",10);
	 	 }

	printf("faulty\n");

	}
rns_to_binary( &Q->x,tempP.x,Bn,pool);
rns_to_binary(&Q->y,tempP.y,Bn,pool);
rns_to_binary(&Q->z,tempP.z,Bn,pool);

ec_point_rns_clear(tempP);
ec_point_rns_clear(R0);
ec_point_rns_clear(R1);
ec_point_clear(V);

for (i=0;i<2*MOD_NUM;i++){
	mpz_clear(val_rns[i]);
	mpz_clear(val_rns_1[i]);
	mpz_clear(mmodp_B[i]);

}
mpz_clear(permut_index);
	gmp_randclear(state);
}

void EC_scalar_mul_no_Base_Perm_no_R_point(ec_point *Q, ec_point P, mpz_t k, rns_base_element_data pool, base base_rand_index[], MmodP modp, mpz_t field_p){

	int Bn[MOD_NUM+1];
	int Bnn[MOD_NUM+1];


	ec_point_rns R0=rns_EC_point_init();
	ec_point_rns R1=rns_EC_point_init();
	//ec_point_rns R2=rns_EC_point_init();



	mpz_t mmodp_B[2*MOD_NUM];
	mpz_t val_rns[2*MOD_NUM];
	mpz_t val_rns_1[2*MOD_NUM];


	int i=0;
	for (i=0; i<2*MOD_NUM;i++){


		mpz_init(mmodp_B[i]);
		mpz_init(val_rns[i]);
		mpz_init(val_rns_1[i]);
	}

	int  permut_index=1;  //choise of RNS bases to be used in all scalar multiplication
										// Currently RNS bases choice 1


	for (i=0;i<MOD_NUM;i++){
		Bn[i]=base_rand_index[permut_index].moduli[i];
		Bnn[i]=base_rand_index[permut_index].m_nn_nums[i];
		}
		Bn[MOD_NUM]=base_rand_index[permut_index].base_index;
		Bnn[MOD_NUM]=69-base_rand_index[permut_index].base_index;

		for (i=0; i<2*MOD_NUM;i++){
			mpz_set_str(val_rns[i],"0",10);
			mpz_set_str(val_rns_1[i],"1",10);
		}
		for (i=0;i<MOD_NUM;i++){
			mpz_set(mmodp_B[Bn[i]],modp.m[Bn[i]]);
			mpz_set(mmodp_B[Bnn[i]],modp.m[Bnn[i]]);
		}


		//base point in rns format

		ec_point V=EC_point_init();

		ec_point_rns tempP=rns_EC_point_init();

			 binary_to_rns(P.x, tempP.x, Bn, pool);
			 binary_to_rns(P.x, tempP.x, Bnn, pool);
			 binary_to_rns(P.y, tempP.y, Bn, pool);
			 binary_to_rns(P.y, tempP.y, Bnn, pool);
			 binary_to_rns(P.z, tempP.z, Bn, pool);
			 binary_to_rns(P.z, tempP.z, Bnn, pool);

			 for (i=0;i<2*MOD_NUM;i++){
			 mpz_set(R1.x[i],tempP.x[i]);
			 mpz_set(R1.y[i],tempP.y[i]);
			 mpz_set(R1.z[i],tempP.z[i]);
			mpz_set_str(R0.x[i],"0",10);
			mpz_set_str(R0.y[i],"1",10);
			mpz_set_str(R0.z[i],"1",10);
		}



//conversion to Montgomery format


		RNS_Montg_mul(R1.x, R1.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.y, R1.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R1.z, R1.z, mmodp_B,Bnn,Bn,pool);

		RNS_Montg_mul(R0.x, R0.x, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.y, R0.y, mmodp_B,Bnn,Bn,pool);
		RNS_Montg_mul(R0.z, R0.z, mmodp_B,Bnn,Bn,pool);


		gmp_printf("----\n");


	for (i=EC_Field_P_bit_lenght-1;i>=0;i--){
		int j;



		if (mpz_tstbit(k,i)==1){

			EdC_point_add(&R0, R1, R0, Bn, Bnn,pool,field_p);

			EdC_point_double(&R1, R1, Bn, Bnn,pool,field_p);

		}
		else {

			EdC_point_add(&R1, R1, R0, Bn, Bnn,pool,field_p);
			EdC_point_double(&R0, R0, Bn, Bnn,pool,field_p);
		}
	}




	RNS_Montg_mul(tempP.x, tempP.x, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.y, tempP.y, mmodp_B,Bnn,Bn,pool);
	RNS_Montg_mul(tempP.z, tempP.z, mmodp_B,Bnn,Bn,pool);

	 EdC_point_add(&tempP, tempP, R0, Bn, Bnn,pool,field_p);

	rns_to_binary(&V.x, R1.x,Bn,pool);
	mpz_sub(V.x,val_rns[0],V.x);
	mpz_mod(V.x,V.x,field_p);
	binary_to_rns(V.x,R1.x,Bn,B_pool);
	binary_to_rns(V.x,R1.x,Bnn,B_pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bn,pool);
	//RNS_sub(R1.x,val_rns,R1.x,Bnn,pool);

	EdC_point_add(&tempP, tempP, R1, Bn, Bnn,pool,field_p);


	if((mpz_cmp(tempP.x[0],val_rns[0])==0) &&
		(mpz_cmp(tempP.x[1],val_rns[1])==0) &&
		(mpz_cmp(tempP.x[2],val_rns[2])==0) &&
		(mpz_cmp(tempP.x[3],val_rns[3])==0) &&
		(mpz_cmp(tempP.x[4],val_rns[4])==0) &&
		(mpz_cmp(tempP.x[5],val_rns[5])==0) &&
		(mpz_cmp(tempP.x[6],val_rns[6])==0) &&
		(mpz_cmp(tempP.x[7],val_rns[7])==0))  //&&

		{




			RNS_Montg_mul(tempP.x, R0.x,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.y, R0.y,val_rns_1,Bn,Bnn,pool);
			RNS_Montg_mul(tempP.z, R0.z,val_rns_1,Bn,Bnn,pool);

	}
	else
	{ for (i=0; i<2*MOD_NUM;i++){
		mpz_set_str(tempP.x[i],"0",10);
		mpz_set_str(tempP.y[i],"1",10);
		mpz_set_str(tempP.z[i],"1",10);
	 	 }

	printf("faulty\n");

	}
rns_to_binary( &Q->x,tempP.x,Bn,pool);
rns_to_binary(&Q->y,tempP.y,Bn,pool);
rns_to_binary(&Q->z,tempP.z,Bn,pool);

ec_point_rns_clear(tempP);
ec_point_rns_clear(R0);
ec_point_rns_clear(R1);
ec_point_clear(V);

for (i=0;i<2*MOD_NUM;i++){
	mpz_clear(val_rns[i]);
	mpz_clear(val_rns_1[i]);
	mpz_clear(mmodp_B[i]);

}

}
