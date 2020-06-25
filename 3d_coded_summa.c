#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mkl.h"


void print_node_msg(MPI_Comm comm, char* msg) { 
	int rank;
	MPI_Comm_rank(comm, &rank); 
	printf("Node %d: %s \n", rank, msg); 
}

void print_matrix_d(double* matrix, int dim1, int dim2){
	int i, j;
	for (i = 0; i < dim1; i++){
		for(j =0; j<dim2; j++){
			printf("%10.3f ", matrix[i*dim2 + j]);
		}
		printf("\n");
	}
	printf("\n");
}

int * rand_perm(int n, int k) {
	int *a = malloc(n*sizeof(int));
	int i;
	for (i = 0; i < n; i++){
		a[i] = i;
	}
	for (i = n-1; i > 0; i--) {
		int j = rand() % (i+1);
		int temp = a[j];
		a[j] = a[i];
		a[i] = temp;
	}

	int *b = malloc(k*sizeof(int));
	for (i=0; i<k; i++) {
		b[i] = a[i];
	}    
	free(a);
	return b;
}



int decode_matdot(int n, int k, int *succ_idx, double *decoding_vec){
	double alpha[n];
	int i, j;

	/* Obatining  Chebyshev nodes */
	for (i=0; i < n; i ++){
		alpha[i] = cos((2*i+1.0)/(2.0*n)*M_PI);
	}                                                         

	/* Vandermonde matrix of the successful indices  */
	double vander[k*k];
	/* Inverse of the Vandermonde matrix */

	for (i=0; i < k; i++){
		for(j=0; j < k; j++){
			vander[i*k+j] = pow(alpha[succ_idx[i]], j);
		}
	}

	/* auxillary variable for LU factoriztion & matrix inversion */
	int ipiv[k];
	/* LU factorizaion of vander*/
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, k, k, vander, k, ipiv );
	/* Inverting the Vandermonde matrix: vander = inv(vander) */
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, k, vander, k, ipiv);

	int idx;
	for (j=0; j < k; j++){
		/* decoding_vec = (int(k/2))-th row of vander */
		idx = (int)(k/2)*k + j;
		decoding_vec[j] = vander[idx];
	}

	return 0;
}




int main(int argc, char **argv) {

	time_t curr_time; 
	time(&curr_time);

	/* For (n,k)x(N,K) coded 2.5D SUMMA */ 
	/* Default values are (4,4)x(4,2) coded 2.5D SUMMA on 24x24 matrix */ 
	int k = 4, n = 4;
	int K = 2, N = 4;
	int m = 48;


	/* Parsing Command line arguments for the parameters n,k,N,K,m */

	int opt;
	bool debug_mode = false;
	bool load_balanced = false;

	while ((opt = getopt(argc, argv, "k:n:K:N:m:dl")) != -1){
		switch (opt) {
			case 'k':
				k = atoi(optarg);
				break;
			case 'n': 
				n = atoi(optarg);
				break;			
			case 'K':
				K = atoi(optarg); 
				break;
			case 'N': 
				N = atoi(optarg);

				break;
			case 'm':
				m = atoi(optarg);
				break;
			case 'd':
				debug_mode = true;
				break;  
			case 'l':
				load_balanced = true;
				break;
			default: 
				exit(-1);
				break;
		}
	}


	if (n < k) {
		printf("n should be greater than or equal to k.\n");
		exit(-1);
	}

	if (N < K) {
		printf("N should be greater than or equal to K.\n");
		exit(-1);
	}

	if (m < (n*N)) { 
		printf("m should be at least n*N.\n");
		exit(-1);
	}

	/***********************
	 ***** INITIALIZATION****
	 ************************/

	MPI_Init(NULL, NULL);


	/* For iterators */
	int i,j, l;

	/* Get the rank and size in the original communicator */
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	double t_begin = 0, t_end = 0; 
	/**  
	 * Split the communicator for x, y, z axes
	 * Define x_comm, y_comm, z_comm and 
	 * x_rank, y_rank, z_rank
	 * 	 	 	 	       **/


	/** This part is hard-coded for a reason.. **/ 


	int p = world_rank;
	int x, y, z, color_x, color_y, color_z;

	/* For (n=8, N=4) [NOTE] npernode=16  */
	/* The logic for n=8 is to do the samething as n=4 case 
	 * for each 4x4 quarter, and add 4 on x or y axes */ 
	if(load_balanced){
		if (8 == n && 4 == N) {
			p = world_rank % 64; 	

			z = p % 4;
			y = (p%16)/4; 
			x = (p-z-4*y)/16;
			x = (x-y-z) % 4; 
			if (x < 0){
				x += 4;
			}

			if (world_rank < 64) {
				;	
			} else if (world_rank < 128) {
				x += 4;	
			} else if (world_rank <192) {
				y += 4;

			} else {
				x += 4;
				y += 4;
			}

			color_x = 4*y+z;
			color_y = 4*x+z;
			color_z = 8*x+y;
		}

		/* For (n=8, N=2) [NOTE] npernode=8 */ 

		else if (8 == n && 2 == N) {
			p = world_rank % 32;

			x = p % 4; 
			z = ((p % 8)-x)/4;
			y = ((p-x-4*z)/8-x-z)%4;
			if (y <0){
				y+= 4; 
			}

			if (world_rank < 32){
				;
			} else if (world_rank < 64){ 
				x += 4;
			} else if (world_rank < 96){
				y += 4;
			} else {
				x += 4;
				y += 4;
			}

			color_x = 8*z+y;
			color_y = 8*z+x;
			color_z = 8*y+x; 
		}


		/* For (n=8, N=6) [NOTE] npernode=24 */	

		else if (8 == n && 6 == N) {
			z = world_rank % N;
			x = (world_rank/N) % 4 + 4*(world_rank/192);
			y = ((world_rank/24) - x - z - 8*(world_rank/192))% n;

			if (y < 0){
				y += n;
			}

			color_x = 8*z+y;
			color_y = 8*z+x;
			color_z = 8*x+y;	
		}


		/* For (n=8, N=8) [NOTE] npernode=16 */	

		else if (8 == n && 8 == N) {
			z = world_rank % N;
			x = (world_rank/N) % 2 + 2*(world_rank/128);
			y = ((world_rank/16)-x-z-8*(world_rank/128))% n;

			if (y < 0){
				y += n;
			}

			color_x = 8*z+y;
			color_y = 8*z+x;
			color_z = 8*x+y;	
		}


		/* For (n=8, N=10) [NOTE] npernode=16 */	

		else if (8 == n && N == 10) {
			y = world_rank % n;
			x = (world_rank/n) % 2 + 2*(world_rank/160);
			z = ((world_rank/16)-x-y-10*(world_rank/160))% N;

			if (z < 0){
				z += N;
			}

			color_x = 8*z+y;
			color_y = 8*z+x;
			color_z = 8*x+y;	
		}


		/* For (n=4, N=4) [NOTE] npernode=16 */	

		else if (4 == n && 4 == N) {
			z = world_rank % n;
			y = ((world_rank-z) % (n*n))/n; 
			x = (world_rank-z-4*y)/16;
			x = (x-y-z) % n; 
			if (x < 0){
				x += n;
			}

			color_x = 4*y+z;
			color_y = 4*z+x;
			color_z = 4*x+y;	
		}

		/* For (n=4, N=2) [NOTE] npernode=8 */

		else if (4 == n && 2 == N) { 
			x = p % 4; 
			z = ((p % 8)-x)/4;
			y = ((p-x-4*z)/8-x-z)%4;
			if (y <0){
				y+= 4; 
			}

			color_x = 4*z+y;
			color_y = 4*z+x;
			color_z = 4*y+x; 
		}

	}	
	/* For a general case */ 
	else {
		color_x = world_rank / n;  
		color_y = world_rank % n + n*(world_rank/(n*n));
		color_z = world_rank % (n*n);  
		x = y = z = world_rank;
	}


	/***************************************************/
	MPI_Comm x_comm, y_comm, z_comm;

	MPI_Comm_split(MPI_COMM_WORLD, color_x, x, &x_comm);
	MPI_Comm_split(MPI_COMM_WORLD, color_y, y, &y_comm);
	MPI_Comm_split(MPI_COMM_WORLD, color_z, z, &z_comm);	

	int x_rank, x_size;
	MPI_Comm_rank(x_comm, &x_rank);	
	MPI_Comm_size(x_comm, &x_size);

	int y_rank, y_size;
	MPI_Comm_rank(y_comm, &y_rank);
	MPI_Comm_size(y_comm, &y_size);  

	int z_rank, z_size;
	MPI_Comm_rank(z_comm, &z_rank);
	MPI_Comm_size(z_comm, &z_size); 
	
	/** 
	 *  Defining Matrix Dimensions 
	 *  Matrices A, B are m-by-m. 
	 *  m should be divisible by n, K, k. 
	 * 	 	 	 	 **/

	t_begin = MPI_Wtime();

	int m_local = m/n;
	int m_MatDot_L0 = m_local*N/K; 
	int m_MatDot = m_local/K; 
	int m_prod = m/k;
	int local_mat_size = m_MatDot * m_prod;

	if(debug_mode) {
		if(world_rank == 0){ 
			printf("m_local: %d, m_MatDot_L0: %d, m_MatDot: %d, m_prod: %d, local_mat_size:%d \n", m_local, m_MatDot_L0, m_MatDot, m_prod, local_mat_size); 
		}
	}

	/** Variables for allocating memory **/
	int memory_idx = 0; 
	double *memory_block = NULL;
	int total_memory = 0; 

	/** Variables for time measurements **/
	double t1 = 0, t2 = 0, t_malloc =0, t_matdot_enc = 0, t_scatter =0, t_prod=0, t_decode=0, t_reduce = 0, t_summa =0; 

	t1 = MPI_Wtime();
	/** Allocating the memory block **/ 
	if (z_rank == 0) { 
		total_memory = 2*m_local*m_local*(1+N/K);
		memory_block = (double *)  mkl_calloc(total_memory, sizeof(double), 64);
	} else {
		total_memory = 4*m_local*m_local*n/k/K + m_prod*m_prod;
		memory_block = (double *) mkl_calloc(total_memory, sizeof(double), 64);
	}

	t2 = MPI_Wtime();
	t_malloc = t2-t1; 	



	double *A_L0 = NULL, *B_L0 = NULL, *A_MatDot_L0 = NULL, *B_MatDot_L0 = NULL;


	if (debug_mode){
		printf("(%d, %d, %d): Assigned %d doubles. \n", x_rank, y_rank, z_rank, total_memory);
	}

	/** Local generation of matrices at Layer 0 **/

	if (z_rank == 0){
		/** Local generation of matrices A_L0 and B_L0 
		 *  (This can be replaced by any other methods to generate matrices) 
		 *  For now, at node P: A_L0[i,j] = (i+j) + 0.01*P,
		 *  B_L0[i,j] = (10*i+j) + 0.01*P
		 * 		 		 		 		 **/

		if (debug_mode){ 
			print_node_msg(MPI_COMM_WORLD, "Starts generating a local matrix"); 
		}

		A_L0 = (double *) &memory_block[memory_idx];
		memory_idx += m_local*m_local; /* Size of matrix A_L0 */
		B_L0 = (double *) &memory_block[memory_idx]; 
		memory_idx += m_local*m_local; /* Size of matrix B_L0 */

		assert(A_L0[0] == 0.0); 
		assert(B_L0[0] == 0.0); 

		for (i =0; i < m_local; i++){
			for (j = 0; j < m_local; j++) {
				A_L0[i*m_local+j] = (double) (i+j) + 0.01*x_rank;
			}
		}

		for (i =0; i < m_local; i++){
			for (j = 0; j < m_local; j++) {
				B_L0[i*m_local+j] = (double) (10.0*i+j) + 0.01*x_rank;
			}
		}

		/** MatDot Encoding at the first layer **/ 
		t1 = MPI_Wtime();

		if (debug_mode){
			print_node_msg(MPI_COMM_WORLD, "Starts MatDot encoding");
		}

		if (N > K) {
			/** We only perform MatDot encoding when N > K 
			 **  N = K means "Uncoded" **/
			A_MatDot_L0 = (double *) &memory_block[memory_idx]; 
			memory_idx += m_MatDot_L0*m_local; /* Size of A_MatDot_L0 */
			B_MatDot_L0 = (double *) &memory_block[memory_idx];
			memory_idx += m_local*m_MatDot_L0; /* At this point, memory should be full. */ 

			assert(0.0 == A_MatDot_L0[0]); 
			assert(0.0 == B_MatDot_L0[0]); 

			assert(0.0 == A_MatDot_L0[0]); 
			assert(0.0 == B_MatDot_L0[0]); 


			double alpha[N]; 
			for (i=0; i < N; i ++){ 
				alpha[i] = cos((2*i+1.0)/(2.0*N)*M_PI);
			}

			/* transposing A_L0 for encoding */ 

			mkl_dimatcopy ('r', 't', m_local, m_local, 1.0, A_L0, m_local, m_local);


			int block_idx;
			for (block_idx = 0; block_idx < N; block_idx++){
				for(l = 0; l < K; l++){ 
					mkl_domatadd ('r' , 'n', 'n', m_MatDot, m_local, 1.0, 
							&A_MatDot_L0[block_idx*m_MatDot*m_local], m_local, pow(alpha[block_idx], l), &A_L0[l*m_MatDot*m_local], 
							m_local, &A_MatDot_L0[block_idx*m_MatDot*m_local], m_local);

					/* B_enc_block = B_enc[block_idx*m_MatDot:(block_idx+1)*m_MatDot, :] */
					mkl_domatadd ('r' , 'n', 'n', m_MatDot, m_local, 1.0, 
							&B_MatDot_L0[block_idx*m_MatDot*m_local], m_local, pow(alpha[block_idx], K-1-l), &B_L0[l*m_MatDot*m_local], 
							m_local, &B_MatDot_L0[block_idx*m_MatDot*m_local], m_local);
				}
			}
		} else {
			A_MatDot_L0 = (double *) &memory_block[memory_idx]; 
			memory_idx += m_local*m_MatDot_L0; /* Size of A_MatDot_L0 */
			B_MatDot_L0 = (double *) &memory_block[memory_idx];
			memory_idx += m_local*m_MatDot_L0; /* At this point, memory should be full. */ 

			mkl_domatcopy ('r', 'n', m_local, m_local, 1.0, A_L0, m_local, A_MatDot_L0, m_local);
			mkl_domatcopy ('r', 'n', m_local, m_local, 1.0, B_L0, m_local, B_MatDot_L0, m_local);

		}

		t2 = MPI_Wtime();
		t_matdot_enc = t2-t1;

	}
	/** Scatter the MatDot-encoded matrices along the z axis
	 *  Each node will receive a column block: A_MatDot (m_local-by-m_MatDot) 
	 *                         a row block: B_MatDot (m_MatDot-by-m_local)
	 * m_MatDot is a block size of one MatDot encoded block: m_MatDot = m_local/K = m_sub_matdot/N. 
	 * Because A_enc is transposed, initially A_MatDot will be m_MatDot-by-m_local,
	 * and it will be transposed back to m_local-by-m_MatDot
	 * 	 	 	 	 	 	 **/


	double *A_MatDot = NULL, *B_MatDot = NULL; 
	/* This will be ignored at Layer 0 later */ 
	/* Just a placekeeper at Layer 0 */ 

	if (z_rank == 0){ 
		memory_idx = 2*m_MatDot * m_prod;
	} else { 
		memory_idx = 0;
	}
	A_MatDot = (double*) &memory_block[memory_idx]; 
	memory_idx += m_local*m_MatDot;
	B_MatDot = (double*) &memory_block[memory_idx]; 
	memory_idx += m_local*m_MatDot;


	assert(NULL != A_MatDot);
	assert(NULL != B_MatDot); 

	if (debug_mode){
		print_node_msg(MPI_COMM_WORLD, "Starts scattering the matrix"); 		
	}

	MPI_Barrier(z_comm);

	t1 = MPI_Wtime();

	MPI_Scatter(A_MatDot_L0, local_mat_size, MPI_DOUBLE, A_MatDot, local_mat_size, MPI_DOUBLE, 0, z_comm); 
	MPI_Scatter(B_MatDot_L0, local_mat_size, MPI_DOUBLE, B_MatDot, local_mat_size, MPI_DOUBLE, 0, z_comm); 
	// if (z_rank == 0) {
	// 	/* Resetting A_MatDot and B_MatDot at the first layer */ 
	// 	/* NOTE: Maybe we can use MPI_IN_PLACE later          */
	// 	A_MatDot = A_MatDot_L0;
	// 	B_MatDot = B_MatDot_L0;
	// }

	/* Transposing A_MatDot back into a column block */
	mkl_dimatcopy('r', 't', m_MatDot, m_local, 1.0, A_MatDot, m_local, m_MatDot);

	t2 = MPI_Wtime(); 
	t_scatter = t2 - t1; 


	/**************************************************/
	/* Completely ignoring the Product encoding steps */
	/* Just assigining memory for it now.             */
	/**************************************************/
	t1 = MPI_Wtime();
	double *A_MatDot_Prod = NULL, *B_MatDot_Prod = NULL; 

	if (z_rank == 0) { 
		memory_idx = 0;
	}	
	A_MatDot_Prod = (double*) &memory_block[memory_idx]; /* dimension: m_MatDot x m_prod */
	memory_idx += m_MatDot * m_prod;
	B_MatDot_Prod = (double*) &memory_block[memory_idx]; 
	memory_idx += m_MatDot * m_prod;

	/* Just copying A_MatDot and B_MatDot into A_MatDot_Prod, B_MatDot_Prod for now */ 
	mkl_domatcopy('r', 'n', m_prod, m_MatDot, 1.0, A_MatDot, m_MatDot, A_MatDot_Prod, m_MatDot);
	mkl_domatcopy('r', 'n', m_MatDot, m_prod, 1.0, B_MatDot, m_prod, B_MatDot_Prod, m_prod);

	t2 = MPI_Wtime();
	t_prod = t2 - t1;

	/**********************************/
	/******* SUMMA computation ********/ 
	/**********************************/
	MPI_Barrier(MPI_COMM_WORLD);	
	t1 = MPI_Wtime();	

	double *A_temp = NULL, *B_temp = NULL, *C = NULL;

	if (z_rank == 0){
		A_temp = (double *) &memory_block[memory_idx];
		memory_idx += local_mat_size;
		B_temp = (double *) &memory_block[memory_idx];
		memory_idx += local_mat_size;
		C = (double *) &memory_block[memory_idx];
		memory_idx += m_prod * m_prod; 
	} else {
		C = (double *)  &memory_block[memory_idx];
		memory_idx = 0; 
		A_temp = (double *) &memory_block[memory_idx];
		memory_idx += local_mat_size;
		B_temp = (double *) &memory_block[memory_idx];
		memory_idx += local_mat_size;
	}


	assert(NULL != C); 

	for (i=0; i < m_prod*m_prod; i++){
		C[i] = 0.0;
	}

	if (debug_mode){
		printf("Beginning SUMMA computation at Node%d = (%d, %d, %d) \n", world_rank, color_x, color_y, color_z); 
		if(world_rank == 13 && m < 100){ 
			print_matrix_d(A_MatDot_Prod, m_prod, m_MatDot); 
			print_matrix_d(B_MatDot_Prod, m_MatDot, m_prod);
		}	
	}

	if(world_rank == 0){
		printf("Starting SUMMA computation. Each node has A,B: (%d x %d), C: (%d x %d)\n", m_prod, m_MatDot, m_prod, m_prod);
	}

	// if (world_rank == 0){
	// 	printf("Start: %d, A_MatDot: %d, B_MatDot:%d \n", (int)(&memory_block[0])/8, (int)A_MatDot/8, (int)B_MatDot/8);
	// 	printf("A_MatDot_Prod: %d, B_MatDot_Prod: %d\n", (int)A_MatDot_Prod/8, (int)B_MatDot_Prod/8);
	// 	printf("IDX: %d, A_temp:%d, B_temp:%d, C: %d, End: %d\n", (int)(&memory_block[memory_idx])/8,
	// 	 (int)A_temp/8, (int)B_temp/8, (int)C/8, (int)(&memory_block[total_memory])/8);
	// }

	int iteration;
	double *A_temp_pointer, *B_temp_pointer;

	for (iteration = 0; iteration < n; iteration++) {
		/* row broacast */

		if (x_rank==iteration) {
			A_temp_pointer = A_temp;
			A_temp = A_MatDot_Prod;
		} 

		assert(NULL != A_temp);

		
		MPI_Barrier(x_comm);
		MPI_Bcast(A_temp, local_mat_size, MPI_DOUBLE, iteration, x_comm);

		if (y_rank==iteration) {
			B_temp_pointer = B_temp;
			B_temp = B_MatDot_Prod; 
		} 

		assert (NULL != B_temp); 

		MPI_Barrier(y_comm);
		MPI_Bcast(B_temp, local_mat_size, MPI_DOUBLE, iteration, y_comm);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
				m_prod, m_prod, m_MatDot, 1.0, A_temp, m_MatDot, B_temp, m_prod, 1.0, C, m_prod);

		if (x_rank == iteration) {
			A_temp  = A_temp_pointer; 
		}

		if (y_rank == iteration) {
			B_temp = B_temp_pointer; 
		}
	}

	t2 = MPI_Wtime(); 
	t_summa = t2-t1;

	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime(); 
	if (N > K) { 
		/** MatDot Decoding 
		 *  Decoding MatDot using MPI Reduce 
		 * 	 	 	 	   **/ 

		/* Example indices of successful workers 
		 ** NOTE: this can change later.. 
		 ** For now, it's assuming that one entire layer has failed! */ 
		/* Recovery threshold of MatDot codes */ 
		int matdot_recov_thrsh = 2*K-1;
		/* Randomly choosing (2K-1) successful layers */
		int *succ_idx = rand_perm(N, matdot_recov_thrsh); 
		/* Retrieving decoding vector for the given successful indices */ 
		double matdot_dec_vec[matdot_recov_thrsh]; 
		decode_matdot(N, matdot_recov_thrsh, succ_idx, matdot_dec_vec); 

		int dec_idx = -1; 
		for (i=0; i < matdot_recov_thrsh; i++){
			if(z_rank == succ_idx[i]) {
				dec_idx = i;
				break; 
			}
		}

		/* Scaling the matrix for decoding 
		 * If the node is on a failed layer, simply scale the matrix by 0. 
		 * If the node is on a successful layer, scale the matrix by appropriate coeffieint
		 * given by matdot_dec_vec.  
		 * 	 	 	 	 	 	 	 	 * */ 
		double scale = 0.0; 
		if (dec_idx >=0) { 
			scale = matdot_dec_vec[dec_idx]; 
		} 

		mkl_dimatcopy('r', 'n', m_prod, m_prod, scale, C, m_prod, m_prod);	

	}

	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	t_decode = t2-t1;


	double *C_decoded = NULL;
	if(0 == z_rank) {
		C_decoded = (double *) &memory_block[memory_idx];	  	
		assert(NULL != C_decoded);
	} 

	t1 = MPI_Wtime();

	MPI_Reduce(C, C_decoded, m_prod*m_prod, MPI_DOUBLE, MPI_SUM, 0, z_comm); 

	t2 = MPI_Wtime();
	t_reduce = t2 - t1;


	if(debug_mode){
		if(0 == world_rank && m < 100) {
			printf("C_decoded at Node %d\n", world_rank);
			print_matrix_d(C_decoded, m_prod, m_prod); 
		}
	}



	t_end = MPI_Wtime(); 

	/** Opening the log file **/
	char filename[50]; 
	sprintf(filename,"logs/n%d_k%d_N%d_K%d_m%d.result", n, k, N, K, m);
	MPI_File fp;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp); 


	/** Writing in the log file **/
	char buf[100]; 

	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	long num_bytes = total_memory*8;

	MPI_Get_processor_name(processor_name, &name_len);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] Total memory: %ld Bytes \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, num_bytes);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] Memory Allocation: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_malloc);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] MatDot Encoding time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_matdot_enc);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] MPI Scatter time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_scatter);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] Product Encoding  time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_prod);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] SUMMA Total time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_summa);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] MatDot Decoding time: %8.3f s \n",
			world_rank, x_rank, y_rank, z_rank, processor_name, t_decode);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);


	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] MPI Reduce time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_reduce);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);

	sprintf(buf, "[Rank %d (%d, %d, %d) at %s] Total Execution  time: %8.3f s \n", 
			world_rank, x_rank, y_rank, z_rank, processor_name, t_end-t_begin);
	MPI_File_write_ordered(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);


	MPI_File_close(&fp);


	MPI_Barrier(MPI_COMM_WORLD);	

	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp); 
	if (world_rank ==0){
		sprintf(buf, "==========================================================================\n");
		MPI_File_write(fp, buf, strlen(buf), MPI_CHAR, MPI_STATUS_IGNORE);
	}

	MPI_File_close(&fp);

	mkl_free(memory_block);

	MPI_Comm_free(&x_comm);
	MPI_Comm_free(&y_comm);
	MPI_Comm_free(&z_comm);

	MPI_Finalize();



	return 0;

}



