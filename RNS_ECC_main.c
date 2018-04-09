/*
 Name        : RNS_ECC.c
 Author      : Apostolos P. Fournaris
 Version     :
 Copyright   : Your copyright notice
 Description : SCA-FA RNS Scalar multiplier
 ============================================================================
 */

/* 
Compile RaspberryPi: gcc RNS.c RNS_ECC.c -o RNS_ECC -lgmp -lpigpio -lrt -lpthread
Compile BBB: gcc RNS.c RNS_ECC.c -o RNS_ECC -lgmp

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include </home/user/workspaceA/RNS_ECC/RNS.h>
#include <sys/time.h>
#include "gmp.h"
#include <termios.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

//#define DEBUG
//#define DEB_MONT

//#include <pigpio.h>

//int main(void) {

int main(int argc, char** argv){

 int fd;
 const char one = '1';
 const char zero = '0';
 //fd = open("/sys/class/gpio/gpio115/value", O_WRONLY);
 //time_t start, tv2;



struct timeval  tv1, tv2;
double start;


	  mpz_t scalar;
	  gmp_randstate_t exp_seed;
	  int n=192;
	  mpz_init(scalar);
//fixed scalar
	if (!strcmp(argv[1],"0")) {
	  mpz_init_set_str (scalar,"6277101735386680763835789423197608153283161448255502417921",10);
	  mpz_out_str(stdout, 10, scalar);
	}
	else {
//random scalar

		  gmp_randinit_default(exp_seed);
		  gettimeofday(&tv1,NULL);
		  gmp_randseed_ui(exp_seed,tv1.tv_usec);
		  mpz_urandomb(scalar, exp_seed, (mp_bitcnt_t ) n);
		  mpz_out_str(stdout, 10, scalar);
		  printf("\n");
	}

	  mpz_t Fp;
	  mpz_init_set_str (Fp,"6277101735386680763835789423207666416083908700390324961279",10);  //p=2^192-2^64-1

	ec_point base_point=EC_point_init();   //base point in projective coordinates


	mpz_set_str(base_point.x, "3645384250405747444443975179707878184422170030891242095163",10);
	mpz_set_str(base_point.y, "2183112785974248518085123806744212598612038813705995087216",10);
	mpz_set_str(base_point.z, "2567337735049316422574991945429690775704826790646489491016",10);



	ec_point random_point=EC_point_init(); //a random point in projective coordinates

		//Random point 1
		//mpz_set_str(random_point.x, "3227322422912483731737114949460790868615472667137656441460",10);
		//mpz_set_str(random_point.y, "1729676664365958710398636854338808642429694153983591117293",10);
		//mpz_set_str(random_point.z, "3445052988600306549398810979482243117851736803280300345792",10);


		//Random point 2
		//mpz_set_str(random_point.x, "1675264768726911020436433596821534541926713021346646324307",10);
  		//mpz_set_str(random_point.y, "4007721353430501891668103935442728489092034015420501022435",10);
  		//mpz_set_str(random_point.z, "1",10);


		//Random point 3
/*
		mpz_set_str(random_point.x,"5756925922632589071825271859744142962297365779024435394019",10);
		mpz_set_str(random_point.y,"293063117661824273195884929415479374553529558314284996144",10);
		mpz_set_str(random_point.z,"5267325044927431309935306240683483612476683036591677935925",10);

*/
		//Double of point 3
		//mpz_set_str(random_point.x,"1038315139677369774954888565161280850897426943238966379223",10);
		//mpz_set_str(random_point.y,"1740598953895877493307093429464081220314659478756423209033",10);
		//mpz_set_str(random_point.z,"3409444143908873623427979825856752702744427443557190588443",10);


	ec_point output_point=EC_point_init(); // the output point in projective coordinates



#ifdef DEBUG
	ec_point_rns base_point_rns=rns_EC_point_init();
	ec_point_rns output_point_rns=rns_EC_point_init();
	ec_point_rns random_point_rns=rns_EC_point_init(); //random point in rns format



		mpz_t MA;
		mpz_init(MA);
#endif

	//rns_base_element_data bl=rns_base_data_init();
//	int permut_index;
	MmodP mp;
	mpz_init(mp.MM);
	int i;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(mp.m[i]);

	}




	generate_base_element_pool(B_pool, Fp, base_index,&mp);
#ifdef DEBUG
	int permut_index=3;
	for (i=0; i<MOD_NUM;i++){
		Bn[i]=base_index[permut_index].moduli[i];
		Bnn[i]=base_index[permut_index].m_nn_nums[i];
	}
	Bn[MOD_NUM]=base_index[permut_index].base_index;
	Bnn[MOD_NUM]=69-base_index[permut_index].base_index;


	binary_to_rns(base_point.x, base_point_rns.x, Bn, B_pool);
	binary_to_rns(base_point.x, base_point_rns.x, Bnn, B_pool);

	binary_to_rns(base_point.y, base_point_rns.y, Bn, B_pool);
	binary_to_rns(base_point.y, base_point_rns.y, Bnn, B_pool);

	binary_to_rns(base_point.z, base_point_rns.z, Bn, B_pool);
	binary_to_rns(base_point.z, base_point_rns.z, Bnn, B_pool);


	binary_to_rns(random_point.x, random_point_rns.x, Bn, B_pool);
	binary_to_rns(random_point.x, random_point_rns.x, Bnn, B_pool);

	binary_to_rns(random_point.y, random_point_rns.y, Bn, B_pool);
	binary_to_rns(random_point.y, random_point_rns.y, Bnn, B_pool);

	binary_to_rns(random_point.z, random_point_rns.z, Bn, B_pool);
	binary_to_rns(random_point.z, random_point_rns.z, Bnn, B_pool);

	mpz_set(MA,B_pool[Bn[0]].modulo);
	for (i=1;i<MOD_NUM;i++){
		mpz_mul(MA,MA, B_pool[Bn[i]].modulo);
	}

	mpz_t val_rns[2*MOD_NUM];
	mpz_t val_rns_1[2*MOD_NUM];

	for (i=0; i<2*MOD_NUM;i++){
			mpz_init(val_rns[i]);
			mpz_init(val_rns_1[i]);
			mpz_set_str(val_rns[i],"0",10);
			mpz_set_str(val_rns_1[i],"1",10);
		}

	mpz_t mmodp_B[2*MOD_NUM];
	//int i=0;
	for (i=0;i<2*MOD_NUM;i++){
		mpz_init(mmodp_B[i]);
	}
	for (i=0;i<MOD_NUM;i++){
		mpz_set(mmodp_B[Bn[i]],mp.m[Bn[i]]);
		mpz_set(mmodp_B[Bnn[i]],mp.m[Bnn[i]]);

	}
	RNS_Montg_mul(base_point_rns.x, base_point_rns.x, mmodp_B,Bnn,Bn,B_pool);
	RNS_Montg_mul(base_point_rns.y, base_point_rns.y, mmodp_B,Bnn,Bn,B_pool);
	RNS_Montg_mul(base_point_rns.z, base_point_rns.z, mmodp_B,Bnn,Bn,B_pool);

	RNS_Montg_mul(random_point_rns.x, random_point_rns.x, mmodp_B,Bnn,Bn,B_pool);
	RNS_Montg_mul(random_point_rns.y, random_point_rns.y, mmodp_B,Bnn,Bn,B_pool);
	RNS_Montg_mul(random_point_rns.z, random_point_rns.z, mmodp_B,Bnn,Bn,B_pool);

	gmp_printf("Base 1 %d, %d, %d, %d  \n", base_index[permut_index].moduli[0],base_index[permut_index].moduli[1], base_index[permut_index].moduli[2], base_index[permut_index].moduli[3]);
	gmp_printf("Base 2 %d, %d, %d, %d  \n", base_index[permut_index].m_nn_nums[0],base_index[permut_index].m_nn_nums[1], base_index[permut_index].m_nn_nums[2], base_index[permut_index].m_nn_nums[3]);


#endif



	//EdC_point_add(&output_point_rns, base_point_rns, random_point_rns, Bn, Bnn, B_pool,Fp);
	//EdC_point_double(&output_point_rns, base_point_rns, Bn, Bnn, B_pool,Fp);
	
	//The following function generates a random point
	Rand_EC_point(&random_point, base_point, B_pool, base_index, mp, Fp);
	//-------------------------

	gmp_printf("Random  Point R= %Zd,  %Zd, %Zd \n -------- \n", random_point.x, random_point.y, random_point.z);

	gettimeofday(&tv1, NULL);

/*
//start triggering for RaspberryPi:
	double start;
	if (gpioInitialise() < 0)
	{
		fprintf(stderr, "pigpio init failed\n");
		return 1;
	}
	
	gpioSetMode(12,PI_OUTPUT);

        	
	start = time_time();
	while ((time_time() - start) < 10000.0)
	{
	gpioWrite(12,1);
	fprintf(stdout, "Triggering\n");
	//time_sleep(0.2);
	EC_scalar_mul(&output_point, base_point, random_point, scalar, B_pool, base_index, mp, Fp);
	//EC_scalar_mul_no_Base_Perm(&output_point, base_point, random_point, scalar, B_pool, base_index, mp, Fp);
	//EC_scalar_mul_no_R_point(&output_point, base_point, scalar, B_pool, base_index, mp, Fp);
	//EC_scalar_mul_no_Base_Perm_no_R_point(&output_point, base_point, scalar, B_pool, base_index, mp, Fp);

	//time_sleep(0.3);
	gpioWrite(12,0);
	//time_sleep(0.3);
//break;
	}
	gpioTerminate();
*/


/* Trigger for BBB: */
//	write(fd, &one, 1);

//	start = tv1.tv_sec;
//	while (start < 10000000000.0)
//	{
	//gpioWrite(12,1);
//	write(fd, &one, 1);
	fprintf(stdout, "Triggering\n");
	//time_sleep(0.2);
	//EC_scalar_mul(&output_point, base_point, random_point, scalar, B_pool, base_index, mp, Fp);
	//EC_scalar_mul_no_Base_Perm(&output_point, base_point, random_point, scalar, B_pool, base_index, mp, Fp);
	//EC_scalar_mul_no_R_point(&output_point, base_point, scalar, B_pool, base_index, mp, Fp);
	EC_scalar_mul_no_Base_Perm_no_R_point(&output_point, base_point, scalar, B_pool, base_index, mp, Fp);


	
//	}
//	write(fd, &zero, 0);

	gettimeofday(&tv2, NULL);

	printf ("Total time = %f seconds \n",
	(double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	(double) (tv2.tv_sec - tv1.tv_sec));
	gmp_printf("Output point L= %Zd \n ------------- \n", output_point.x, output_point.y, output_point.z);

//new code with and without randomization

//	EC_scalar_mul_no_Base_Perm_no_R_point(&output_point, base_point, scalar, B_pool, base_index, mp, Fp);



#ifdef DEBUG


RNS_Montg_mul(output_point_rns.x, output_point_rns.x, val_rns_1,Bn,Bnn, B_pool);
RNS_Montg_mul(output_point_rns.y, output_point_rns.y, val_rns_1,Bn,Bnn, B_pool);
RNS_Montg_mul(output_point_rns.z, output_point_rns.z, val_rns_1,Bn,Bnn, B_pool);

rns_to_binary(&output_point.x, output_point_rns.x,Bn,B_pool);
rns_to_binary(&output_point.y, output_point_rns.y,Bn,B_pool);
rns_to_binary(&output_point.z, output_point_rns.z,Bn,B_pool);
#endif

gmp_printf("Output  Point L= %Zd,  %Zd, %Zd \n -------- \n", output_point.x, output_point.y, output_point.z);
gmp_printf("scalar= %Zd \n",scalar);

//clear(scalar);



	return EXIT_SUCCESS;
}


