
// Utilities and system includes
//#include <shrUtils.h>
//#include "cutil_inline.h"
#include "TNT.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
using namespace std;

// includes, kernels
#include <TNT_kernel.cu>


////////////////////////////////////////////////////////////////////////////////
// declaration, forward
////////////////////////////////////////////////////////////////////////////////
extern "C" //prints messages
void msgs(int);

extern "C" //get file location of DB1
char *getDB1FileLoc();

extern "C" //retrieve DB1 sequences from file
void getDB1( int*, int*, int*, int*, char**, int**, int**, int**, char*, int*, int*);    

extern "C" //print DB1 data to a file
void printD( int*, int*, int*, int*, char**, int**, int**, int**, char*, int*, int*); 

extern "C" //retrieve V sequences from V files
void getV(char*, int*, int*, int*);

extern "C" //retrieve J sequences from J files
void getJ(char*, int*, int*, int*);

extern "C" //print V sequences to a file, to check correctness
void printV(char*, int*, int*);

extern "C" //print J sequences to a file, to check correctness
void printJ(char*, int*, int*);

extern "C" //Get InVivo sequence from file
void getInVivo( int*, int*, unsigned char**); 

extern "C" //Get InVivo sequence from file
void printInVivo( int*, int*, unsigned char**); 

extern "C" //Get number of pairs in each VJ combinations
void getNum_VJ_Pairs( int, int, int*, int*); 

int	//determine number of threads per block to use
threads_per_block(int);

void //reduce the result on host machine
hostReduce(int, int, int, int, unsigned int*);

int //determine which kernel we need to run
whichKernel(int);

void //*$*$*$*$gets starting and ending V J indices for a given pair. Not currently set up for palandromic sequences*$*$*$*$
getVJ_Indices(int, int, int*, int*, int*, int*);

extern "C" //print the InVivo Results
void print_InVIvo_Results(int, int, int, int, unsigned int*, int);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	int i, j;
    clock_t t1, t2, t3, t4;
    float diff1, diff2;
	unsigned int trans_Size;	//size of memory to transfer from GPU to host for intermediate results
    int threads;						//number of threads per thread-block
    int blocks;							//number of thread-blocks per grid
	unsigned int total_threads;			//total number of threads in grid
	int n_val;							//n value within a loop
	int h_V_Begin, h_V_End, h_J_Begin, h_J_End; //used for indexing in kernel

	////////////////////////////////////////////////////////////////////////////
	//Global arrays for chewback sequences
	char h_chewArrV[45000];
	char h_chewArrJ[45000];
	size_t size = 45000 * sizeof(char);
	int h_counterIndex = 0;
	////////////////////////////////////////////////////////////////////////////

	//print error message if n is out of bounds. Function msgs(int) found in TNT_gold.cpp
	if( N > 12 || N < 0){msgs(0); exit(1);} 
	
	//Use device with highest Gflops/s performance
	//cudaSetDevice(0);// cutGetMaxGflopsDeviceId() );
    //cudaSetDevice(cutGetMaxGflopsDeviceId());

	//Get the location of the DB1 file. Function found in TNT_gold.cpp
	gene_loc_cp = getDB1FileLoc(); //printf("%s\n", gene_loc_cp);

	//Retrieve data for DB1 from file. Function found in TNT_gold.cpp
	getDB1( &numDB1Unique, 
		    &DB1size, 
			&memSizeDB1, 
			&memSizeUniqueDB1, 
			&h_DB1_cp, 
		    &h_numOccurrenceDB1_ip, 
			&h_numUniqueCharDB1_ip, 
			&h_DB1_base_ip, 
			gene_loc_cp, 
			&h_numCharFullDB1, 
			&h_D1Occur);

	//print out all D sequences to test file. Function found in TNT_gold.cpp
	printD( &numDB1Unique, 
		    &DB1size, 
			&memSizeDB1, 
			&memSizeUniqueDB1, 
			&h_DB1_cp, 
		    &h_numOccurrenceDB1_ip, 
			&h_numUniqueCharDB1_ip, 
			&h_DB1_base_ip, 
			gene_loc_cp, 
			&h_numCharFullDB1, 
			&h_D1Occur);

	//Retrieve data for V from files. Function found in TNT_gold.cpp
	getV(h_V_cp, 
		 h_numUniqueCharV_ip, 
		 h_V_base_ip, 
		 numVUnique_ip);

	//print out all V sequences to test file. Function found in TNT_gold.cpp
	printV(h_V_cp, 
		   h_numUniqueCharV_ip,
		   h_V_base_ip);

	//Retrieve data for J from files. Function found in TNT_gold.cpp
	getJ(h_J_cp, 
		 h_numUniqueCharJ_ip, 
		 h_J_base_ip, 
		 numJUnique_ip);

	//print out all J sequences to test file. Function found in TNT_gold.cpp
	printJ(h_J_cp, 
		   h_numUniqueCharJ_ip,
		   h_J_base_ip);

	//Retrieve the InVivo data
	getInVivo(	&InVivo_memSize64,
				&h_num_InVivo, 
				&h_InVivo_cp64);

	//Prints the InVivo data to a file
	printInVivo(&InVivo_memSize64,
				&h_num_InVivo, 
				&h_InVivo_cp64);

	getNum_VJ_Pairs( NUM_V_FILES, 
					 NUM_J_FILES, 
					 VJ_Pairs_ip,
					 &VJ_Largest); 



	//store number of VJ pair information in constant memory
	cudaMemcpyToSymbol(const_d_VJ_Pairs, VJ_Pairs_ip, NUM_V_FILES*NUM_J_FILES*sizeof(int));

	//find and store base address of each VJ pair file in constant memory
	VJ_Pair_Base_ip[0] = 0;
	j = (NUM_V_FILES*NUM_J_FILES);
	for(i=1;i<j; i++){
		VJ_Pair_Base_ip[i] = VJ_Pairs_ip[i-1] + VJ_Pair_Base_ip[i-1];
	}
	cudaMemcpyToSymbol(const_VJ_Pair_Base, VJ_Pair_Base_ip, NUM_V_FILES*NUM_J_FILES*sizeof(int));


	//Store DB1, V, and J Sequences to constant memory
	//DB1 Sequences
	cudaMemcpyToSymbol(const_d_DB1, h_DB1_cp, memSizeDB1_1);
	cudaMemcpyToSymbol(const_d_DB1_base, h_DB1_base_ip, memSizeDB1_2);
	cudaMemcpyToSymbol(const_d_numOccurrenceDB1, h_numOccurrenceDB1_ip, memSizeDB1_2);
	cudaMemcpyToSymbol(const_d_numUniqueCharDB1, h_numUniqueCharDB1_ip, memSizeDB1_2);
	//V Sequences
	cudaMemcpyToSymbol(const_d_V, h_V_cp, memSizeV);
	cudaMemcpyToSymbol(const_d_V_base, h_V_base_ip, memSizeOffsetV);
	cudaMemcpyToSymbol(const_d_numUniqueCharV, h_numUniqueCharV_ip, memSizeOffsetV);
	//J Sequences
	cudaMemcpyToSymbol(const_d_J, h_J_cp, memSizeJ);
	cudaMemcpyToSymbol(const_d_J_base, h_J_base_ip, memSizeOffsetJ);
	cudaMemcpyToSymbol(const_d_numUniqueCharJ, h_numUniqueCharJ_ip, memSizeOffsetJ);

	//number of times a full chewback of DB1 or DB2 occurs
	cudaMemcpyToSymbol(c_DB_Full_Chew_Occur, &h_D1Occur, sizeof(int));

	//////////////////////////////////////////////////////////////////////////////////////////////
	//Dynamically allocate device chewback arrays and row iterator
	char *d_chewArrV;
	char *d_chewArrJ;
	cudaMalloc((void **)&d_chewArrV, size);
	cudaMalloc((void **)&d_chewArrJ, size);
	cudaMemset(&d_counterIndex, 0, sizeof(int));
	cudaMemcpy(d_chewArrV, h_chewArrV, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_chewArrJ, h_chewArrJ, size, cudaMemcpyHostToDevice);
	/////////////////////////////////////////////////////////////////////////////////////////////

	//allocate memory on GPU for InVivo Sequences
	char* d_InVivo_cp64;  
	//cutilSafeCall( cudaMalloc((void**) &d_InVivo_cp64,  InVivo_memSize64) );
	cudaMalloc((void**) &d_InVivo_cp64,  InVivo_memSize64) ;
	
	//transfer InVivo data to the GPU
	//cutilSafeCall(cudaMemcpy(d_InVivo_cp64,  h_InVivo_cp64,  InVivo_memSize64,  cudaMemcpyHostToDevice) );
	cudaMemcpy(d_InVivo_cp64,  h_InVivo_cp64,  InVivo_memSize64,  cudaMemcpyHostToDevice) ;
	printf("\nInVivo64 Memsize = %d\n\n", InVivo_memSize64);
        fflush(stdout);

	printf("Largest number of VJ pairs %d\n", VJ_Largest);
        fflush(stdout);
	//getting ready to allocate memory for results on host and the GPU
	unsigned int numCombN = (unsigned int)pow(4.0, N); //number of unique n combinations, also equal to total number of threads. Assume n 12 is largest

	//Determine largest possible memory size we will need
	unsigned int result_Memsize;
	if(N == 12)
		result_Memsize = (numCombN/512) * VJ_Largest * sizeof(unsigned int); //total threads * num results per thread * size of int
	else if( (N > 5) && (N < 12) )
		result_Memsize = (numCombN/64) * VJ_Largest * sizeof(unsigned int); 
	else if ( (N > 2) && (N < 6) )
		result_Memsize = (numCombN/64)  * VJ_Largest * sizeof(unsigned int); 
	else
		result_Memsize = VJ_Largest * sizeof(unsigned int); //only one thread block, with one result per sequence in this case

	printf("Memory needed on GPU for iteration results: %1.2f MB\n\n", (float) result_Memsize/1048576 );
        fflush(stdout);

	//allocate memory on the GPU and Host to hold largest set of results. Each thread will store 
	//it's own results cpu will perform reduction while GPU creates new results
	unsigned int *d_Results;
	//cutilSafeCall( cudaMalloc((void**) &d_Results, result_Memsize) );
	cudaMalloc((void**) &d_Results, result_Memsize);
	
	//allocate memory for results for n of largest size on host
	unsigned int *h_Results = (unsigned int *) malloc(result_Memsize);
    t1 = clock();
	for(n_val = N_BEGIN; n_val < N+1; n_val++){	    //iterate through each N value
        t3 = clock();
		for(i=0; i<NUM_V_FILES; i++){				//iterate through each of the VJ pair file combinations
			for(j=0; j<NUM_J_FILES; j++){	
                
				getVJ_Indices(i, j, &h_V_Begin, &h_V_End, &h_J_Begin, &h_J_End);

				cudaMemcpyToSymbol(c_V_Begin, &h_V_Begin, sizeof(int));
				cudaMemcpyToSymbol(c_J_Begin, &h_J_Begin, sizeof(int));
				cudaMemcpyToSymbol(c_V_End, &h_V_End, sizeof(int));
				cudaMemcpyToSymbol(c_J_End, &h_J_End, sizeof(int));

				cudaMemcpyToSymbol(c_n, &n_val, sizeof(int));	//update which n value we are working with in constant memory
				cudaMemcpyToSymbol(c_Vnum, &i,  sizeof(int));	//update current V File in constant memory
				cudaMemcpyToSymbol(c_Jnum, &j,  sizeof(int));	//update current J File in constant memory

				total_threads = (unsigned int)pow(4.0, n_val);	//number of unique n combinations
				threads = threads_per_block(n_val);				//determine number of threads per block
				blocks = total_threads/threads;

				dim3 thread(threads, 1, 1);						//telling GPU threads per block
				dim3 grid(blocks, 1, 1);						//telling GPU total number of blocks

				//----------------------------------//
				//-          Execute Kernel        -//
				//----------------------------------//
				TNT_kernel_InVivo64 <<< grid, thread >>>(d_Results, d_InVivo_cp64, d_chewArrV, d_chewArrJ);
				cudaDeviceSynchronize();

				//how much memory do we need to transfer for the result. i*NUM_J_FILES + j; this is the current file out of 240 we are on
				trans_Size = VJ_Pairs_ip[i * NUM_J_FILES + j] * blocks * sizeof(unsigned int);

				//transfer new set of results to host machine
//				cutilSafeCall(cudaMemcpy(h_Results, d_Results, trans_Size, cudaMemcpyDeviceToHost) );
				cudaMemcpy(h_Results, d_Results, trans_Size, cudaMemcpyDeviceToHost) ;
				
				////////////////////////////////////////////////////////////////////////////////////
				//Transfer memory from device back to host for chewback array
				cudaMemcpy(h_chewArrV, d_chewArrV, size, cudaMemcpyDeviceToHost);
				cudaMemcpy(h_chewArrJ, d_chewArrJ, size, cudaMemcpyDeviceToHost);
				cudaMemcpyFromSymbol(&h_counterIndex, d_counterIndex, sizeof(int), 0, cudaMemcpyDeviceToHost);
				////////////////////////////////////////////////////////////////////////////////////

				//reduce data from GPU result array and store on host result array. 
				hostReduce(n_val, blocks, threads, VJ_Pairs_ip[i * NUM_J_FILES + j], h_Results);
				print_InVIvo_Results(n_val, i, j, trans_Size/blocks/sizeof(int), h_Results, VJ_Pair_Base_ip[i*NUM_J_FILES + j]);
				
			}
		printf("Program %1.2f percent complete for n = %d\n", (float)(i+1)/NUM_V_FILES * 100, n_val);
                fflush(stdout);
		}
        t4 = clock();
        diff2 = (float)t4 - (float)t3;
        cout << "Total time for N = " << n_val << " is = " << diff2/CLOCKS_PER_SEC << " seconds\n" << std::flush;
	}
    t2 = clock();
    diff1 = (float)t2 - (float)t1;
    cout << "Total time for program = " << diff1/CLOCKS_PER_SEC << " seconds\n" << std::flush;

	//int k;
	//Print Chewed back V Sequences
	//cout << "Printed V chewback sequences:" << endl;
	//for (int j=h_V_Begin; j<h_V_End; j++) {
	//	k = h_V_base_ip[j];
	//	for (int i=0; i<h_numUniqueCharV_ip[j]; i++){
	//		cout << h_V_cp[i] << std::flush;
	//		k++;
	//	}
	//	cout << endl;
	//}
	//Print Chewed back J Sequences
	//cout << "Printed J chewback sequences:" << endl;
	//for (int j=h_J_Begin; j<h_J_End; j++) {
	//	k = h_J_base_ip[j];
	//	for (int i=0; i<h_numUniqueCharJ_ip[j]; i++){
	//		cout << h_J_cp[i] << std::flush;
	//		k++;
	//	}
	//	cout << endl;
	//}

	//Print Chewed back V Sequences
	cout << "Printed V chewback sequences:" << endl;
	for (int i = 0; i < 1500; i++) {
		for (int j = 0; j < 30; j++) {
			cout << h_chewArrV[i*30+j] << std::flush;
		}
		cout << endl;
	}
	//Print Chewed back J Sequences
	cout << "Printed J chewback sequences:" << endl;
	for (int i = 0; i < 1500; i++) {
		for (int j = 0; j < 30; j++) {
			cout << h_chewArrJ[i*30+j] << std::flush;
		}
		cout << endl;
	}

	//Free up the GPU memory
    //cutilSafeCallNoSync(cudaFree(d_Results));
	//cutilSafeCallNoSync(cudaFree(d_InVivo_cp64));
	cudaFree(d_Results);
	cudaFree(d_InVivo_cp64);
	cudaFree(d_chewArrV);
	cudaFree(d_chewArrJ);

	//Free up host memory
	free(h_Results);
	free(h_DB1_cp);
	free(h_numOccurrenceDB1_ip);
	free(h_numUniqueCharDB1_ip);
	free(h_DB1_base_ip);
	free(h_V_cp);
	free(h_numUniqueCharV_ip);
	free(h_V_base_ip);
	free(numVUnique_ip);
	free(h_J_cp);
	free(h_numUniqueCharJ_ip);
	free(h_J_base_ip);
	free(numJUnique_ip);
	free(h_InVivo_cp64);
	free(VJ_Pairs_ip);
	free(VJ_Pair_Base_ip);
	//free(h_chewArrV);
	//free(h_chewArrJ);


//    shrEXIT(argc, (const char**)argv);
}



////////////////////////////////////////////////////////////
//Reduce our result on host machine
////////////////////////////////////////////////////////////
void 
hostReduce(int n, int blocks, int threads, int num_elements, unsigned int* h_Results){

	int i, j;
	int sum;

	for(i = 0; i < num_elements; i++){
		sum = 0;
		for(j = 0; j < blocks; j++){
			sum += h_Results[i*blocks + j];
		}
		h_Results[i] = sum;
	}
}

////////////////////////////////////////////////////////////
//determine how many threads we should use per thread-block
////////////////////////////////////////////////////////////
int
threads_per_block(int n_val){

	int threads;

	//determine how many threads we should use
	switch(n_val){
		   case 0:threads =    1; break;
		   case 1:threads =    4; break; 
		   case 2:threads =   16; break;
		   case 3:threads =   64; break;
		   case 4:threads =   64; break;
		   case 5:threads =   64; break;
		   case 6:threads =  128; break;
		   case 7:threads =  64; break; 
		   case 8:threads =  128; break; //128
		   case 9:threads =  128; break; //128
		   case 10:threads = 128; break; //128
		   case 11:threads = 128; break;
		   case 12:threads = 512; break;
	}

	return threads;
}


void
getVJ_Indices(int V_File, int J_File, int* h_V_Begin, int* h_V_End, int* h_J_Begin, int* h_J_End){
	
	switch(V_File){
		case 0:  *h_V_Begin =   0; *h_V_End =  18; break; //h_V_End = last sequence in a V file + 1
		case 1:  *h_V_Begin =  18; *h_V_End =  36; break;
		case 2:  *h_V_Begin =  36; *h_V_End =  54; break;
		case 3:  *h_V_Begin =  54; *h_V_End =  72; break;
		case 4:  *h_V_Begin =  72; *h_V_End =  90; break;
		case 5:  *h_V_Begin =  90; *h_V_End = 106; break;
		case 6:  *h_V_Begin = 106; *h_V_End = 122; break;
		case 7:  *h_V_Begin = 122; *h_V_End = 139; break;
		case 8:  *h_V_Begin = 139; *h_V_End = 156; break;
		case 9:  *h_V_Begin = 156; *h_V_End = 173; break;
		case 10: *h_V_Begin = 173; *h_V_End = 190; break;
		case 11: *h_V_Begin = 190; *h_V_End = 208; break;
		case 12: *h_V_Begin = 208; *h_V_End = 226; break;
		case 13: *h_V_Begin = 226; *h_V_End = 244; break;
		case 14: *h_V_Begin = 244; *h_V_End = 261; break;
		case 15: *h_V_Begin = 261; *h_V_End = 276; break;
		case 16: *h_V_Begin = 276; *h_V_End = 294; break;
		case 17: *h_V_Begin = 294; *h_V_End = 312; break;
		case 18: *h_V_Begin = 312; *h_V_End = 327; break;
		case 19: *h_V_Begin = 327; *h_V_End = 342; break;
	}

	switch(J_File){
		case 0:  *h_J_Begin =   0; *h_J_End =  21; break; //h_J_End = last sequence in a J file + 1 21
		case 1:  *h_J_Begin =  21; *h_J_End =  42; break; //+21
		case 2:  *h_J_Begin =  42; *h_J_End =  65; break; //+23
		case 3:  *h_J_Begin =  65; *h_J_End =  89; break; //+24
		case 4:  *h_J_Begin =  89; *h_J_End = 112; break; //+23
		case 5:  *h_J_Begin = 112; *h_J_End = 138; break; //+26
		case 6:  *h_J_Begin = 138; *h_J_End = 161; break; //+23
		case 7:  *h_J_Begin = 161; *h_J_End = 185; break; //+24+
		case 8:  *h_J_Begin = 185; *h_J_End = 207; break; //22
		case 9:  *h_J_Begin = 207; *h_J_End = 229; break; //+22
		case 10: *h_J_Begin = 229; *h_J_End = 251; break; //+22
		case 11: *h_J_Begin = 251; *h_J_End = 271; break; //+20
	}

}


