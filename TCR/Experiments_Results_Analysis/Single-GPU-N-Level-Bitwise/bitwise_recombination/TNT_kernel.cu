#ifndef _TNT_KERNEL_H_
#define _TNT_KERNEL_H_

#include <stdio.h>

//DB information
__constant__ int const_numDB1 = 169;				//total number of DB1 chewbacks
__constant__ char const_d_DB1[425];				//constant memory allocation for DB chewbacks minus full chew back
//__constant__ char const_d_DB1[1448];				//constant memory allocation for DB chewbacks minus full chew back
__constant__ int const_d_DB1_base[169];				//constant memory contains location of each starting sequence in d_DB1 
__constant__ int const_d_numOccurrenceDB1[169];		//Number of ways for particular DB1 chewback 
__constant__ int const_d_numUniqueCharDB1[169];		//number of characters in a unique occurence of DB1
//V information
__constant__ int const_numV = 342;					//total number of V sequences in all V files
__constant__ char const_d_V[913];					//holds all V chewback sequences
//__constant__ char const_d_V[3107];					//holds all V chewback sequences
__constant__ int const_d_V_base[342];				//contains the starting index of each V sequence
__constant__ int const_d_numUniqueCharV[342];		//number of characters in a unique occurence of V
//J information
__constant__ int const_numJ = 271;					//total number of J sequences in all J files
__constant__ char const_d_J[908];					//holds all J chewback sequences
//__constant__ char const_d_J[3210];					//holds all J chewback sequences
__constant__ int const_d_J_base[271];				//contains the starting index of each J sequence
__constant__ int const_d_numUniqueCharJ[271];		//number of characters in a unique occurence of J


__constant__ int c_DB_Full_Chew_Occur;				//current V sequence
__constant__ int c_Vnum;							//current V sequence
__constant__ int c_Dnum;							//current D sequence
__constant__ int c_Jnum;							//current J sequence
__constant__ int c_n;								//current n value

__constant__ int c_V_Begin;							//Base index for V sequences
__constant__ int c_V_End;							//End index for V sequences
__constant__ int c_J_Begin;							//Base index for J sequences
__constant__ int c_J_End;							//End index for J sequences

__constant__ int const_d_VJ_Pairs[NUM_V_FILES*NUM_J_FILES];
__constant__ int const_VJ_Pair_Base[NUM_V_FILES*NUM_J_FILES];
__constant__ int c_NUM_V_FILES = 20;
__constant__ int c_NUM_J_FILES = 12;

/////////////////////////////////////////////////
//kernel for 64 threads or less
/////////////////////////////////////////////////
__global__ void
TNT_kernel_InVivo64(unsigned int* d_Results, char* d_InVivo_cp64)	
{

	volatile __shared__ char iterSeq_sm[16]; //the thread block size we will use for this kernel is 64
	volatile __shared__ int result_sm[128];  //the max thread-block size

	//The four possible bases
//	char base[4] = {'A', 'T', 'G', 'C'};
	unsigned char base[4] = { 0x00, 0x01, 0x02, 0x03 };	
	__shared__ int arraylocal[9*128];//{Vnum,Jnum, sum, num_Seqs,gl_index, pairBase, seqLen, g_tid, length};
//	char nSeq[12];					//will hold a single n combination
	 arraylocal[0*128 + threadIdx.x]= c_Vnum;				//current V file --Vnum 
	 arraylocal[1*128 + threadIdx.x] = c_Jnum;				//current J file --Jnum
	int sh_index;					//used as a shared memory index
	 //sum;						//holds an iterative sum for result
	
	//obtain a unique global index for each thread in the grid
	arraylocal[7*128 + threadIdx.x] = blockIdx.x * blockDim.x + threadIdx.x; // unsigned int  g_tid
/*
	if(c_n){
		nSeq[0]  = base[g_tid%4];					//n = 1
		nSeq[1]  = base[(g_tid+(g_tid/4))%4];		//n = 2
		nSeq[2]  = base[(g_tid+(g_tid/16))%4];		//n = 3
		nSeq[3]  = base[(g_tid+(g_tid/64))%4];		//n = 4
		nSeq[4]  = base[(g_tid+(g_tid/256))%4];		//n = 5
		nSeq[5]  = base[(g_tid+(g_tid/1024))%4];	//n = 6
		nSeq[6]  = base[(g_tid+(g_tid/4096))%4];	//n = 7
		nSeq[7]  = base[(g_tid+(g_tid/16384))%4];	//n = 8
		nSeq[8]  = base[(g_tid+(g_tid/65536))%4];	//n = 9
		nSeq[9]  = base[(g_tid+(g_tid/262144))%4];	//n = 10
		nSeq[10] = base[(g_tid+(g_tid/1048576))%4]; //n = 11
		nSeq[11] = base[(g_tid+(g_tid/4194304))%4];	//n = 12
	}
*/

	unsigned char encNSeq[3];
	if(c_n){
		encNSeq[0] = ( base[ arraylocal[7*128 + threadIdx.x] % 4 ] << 6 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 4 )) % 4 ] << 4 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 16 )) % 4 ] << 2 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 64 )) % 4] );
		encNSeq[1] = ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] /   256 )) % 4 ] << 6 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] /   1024 )) % 4 ] << 4 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] /    4096 )) % 4 ] << 2 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] /   16384 )) % 4 ] );
		encNSeq[2] = ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 65536 )) % 4 ] << 6 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 262144 )) % 4 ] << 4 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 1048576 )) % 4 ] << 2 ) | ( base[( arraylocal[7*128 + threadIdx.x] + ( arraylocal[7*128 + threadIdx.x] / 4194304 )) % 4 ] );
	}	

	//get the number of InVivo VJ sequences we need to go through
	arraylocal[3*128 + threadIdx.x] = const_d_VJ_Pairs[arraylocal[0*128 + threadIdx.x]*12 + arraylocal[1*128 + threadIdx.x]];				//multiply by 12. Number of J files. --num_Seqs

	//int whichSeq; //which sequence is our current thread-block working on in the scope of current VJ
	//int seqLen;	  //length of our current sequence
	arraylocal[5*128 + threadIdx.x] = const_VJ_Pair_Base[arraylocal[0*128 + threadIdx.x]*12 + arraylocal[1*128 + threadIdx.x]] * 16; //The base address for a given VJ pair --pairBase

	//iterate through all InVivo combinations for current VJ pair
	for( int i = 0; i < arraylocal[3*128 + threadIdx.x]; i++ ) {

		result_sm[ threadIdx.x ] = 0;
		arraylocal[2*128 + threadIdx.x] = 0; //reset our result

		__syncthreads();

		//store an InVivo combination into the shared memory "iterResults_sm[]"
		if( blockDim.x < 16 ) {										//iter through VJ seq if block dim < 64. There's only 1 block
			for( int j = 0; j < ( 16 / blockDim.x ); j++ ) {		//iterations = sequence allocation / block size
				//int k = j * blockDim.x + threadIdx.x;				//create local SM index
				arraylocal[4*128 + threadIdx.x] = ( arraylocal[5*128 + threadIdx.x]  + i * 16 ) + j * blockDim.x + threadIdx.x;			//create global memory index
				iterSeq_sm[ j * blockDim.x + threadIdx.x ] = d_InVivo_cp64[ arraylocal[4*128 + threadIdx.x] ];    	//read the current InVivo sequence from the global memory
			}
		}
		else{	//only threads < 16 will read inVivo data
			if( threadIdx.x < 16 ) {
				arraylocal[4*128 + threadIdx.x] = ( arraylocal[5*128 + threadIdx.x]  + i * 16 ) + threadIdx.x;		//create global memory index
				iterSeq_sm[ threadIdx.x ] = d_InVivo_cp64[ arraylocal[4*128 + threadIdx.x] ];		//read the current InVivo sequence from the global memory
			}
		}

		//if(blockDim.x > 1) 
		__syncthreads();

		//get the length of current sequence for all threads in current thread-block
		// iterSeq_sm[0] = [7:3] V, [2:0] MS-3bits J
		// iterSeq_sm[1] = [7] LS-1bit J, [6:0] Length of InVivo Sequence
		// iterSeq_sm[2] to iterSeq_sm[15] = InVivo code

		unsigned char tempChar0, tempChar1, tempChar2;

		char getChar;
		getChar = iterSeq_sm[1] & 0x7F;				// Mask [7]th bit, as Length = [6:0] ( 0b01111111 = 0x7F )
		arraylocal[6*128 + threadIdx.x] = (int)getChar;						// seqLen = Length of InVivo sequence --seqLen

//		printf("i: %d, seqLen: %d, iterSeq_sm: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", i, seqLen, iterSeq_sm[0], iterSeq_sm[1], iterSeq_sm[2], iterSeq_sm[3], iterSeq_sm[4], iterSeq_sm[5], iterSeq_sm[6], iterSeq_sm[7], iterSeq_sm[8], iterSeq_sm[9], iterSeq_sm[10], iterSeq_sm[11], iterSeq_sm[12], iterSeq_sm[13], iterSeq_sm[14], iterSeq_sm[15]);
//		printf("c_V_Begin: %d, c_V_End: %d \n", c_V_Begin, c_V_End );

		//set our shared memory index to the base of the sequence characters in shared array
		int byteCount = 0;
		int binCount = 0;
		int accuBinCount = 0;
		int k;
		
		int tempAccuBinCount = 0;
		int tempShIndex = 0;

		//int length;				// Length of each generated sequence
		bool Vmatch = true;		// Is there a V sequence match?
		bool seqMatch = true;		// Is the entire sequence a match?

		//////////////////////////////////////////////////////////////////////////////////
		//First compare our InVivo Sequences containing 'V' = Vn, VnJ, VnDn, VnDnJ
		//////////////////////////////////////////////////////////////////////////////////

		// Get all the V sequences:
		for ( int Vindx = c_V_Begin; Vindx < c_V_End; Vindx++ ) {
			
			// Initialize values:
			Vmatch = true;
			seqMatch = true;

			/////////////////////////////////////////////////////////
			//Compare InVivo Sequence to Vn comb with D and J chewed 
			/////////////////////////////////////////////////////////

			arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharV[ Vindx ] + c_n; //length
			if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] ) {
//				printf("Length matches: i: %d, Vindx: %d, c_n: %d, Vlen: %d, length: %d, seqLen: %d\n", i, Vindx, c_n, const_d_numUniqueCharV[ Vindx ], length, seqLen );
				
				sh_index = 2;													// InVivo sequence starts from iterSeq_sm[2]
				accuBinCount = 0;

				byteCount = (( const_d_numUniqueCharV[ Vindx ] ) / 4 );			// Calculates the full bytes of V
				binCount = ( const_d_numUniqueCharV[ Vindx ] % 4 ) * 2;			// Calculates the overflow bits of V
				k = const_d_V_base[ Vindx ];									// Starting address of V sequence

				// Compare the full bytes of V with InVivo sequence
				for ( int m = 0; m < byteCount; m++ ) {
					tempChar0 = const_d_V[ k ];
					tempChar1 = iterSeq_sm[ sh_index ];
					if ( tempChar0 != tempChar1 ) {
						Vmatch = false;
						break;
					}
					sh_index++;
					k++;
				}

				// V does not match => get next V sequence => go to line 158 ( V for loop )
				if ( Vmatch == false ) {
					continue;
				}

				// Compare the overflow bits of V with Invivo sequence
				tempChar0 = ((( const_d_V[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
				tempChar1 = ((( iterSeq_sm[ sh_index ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
				if ( tempChar0 != tempChar1 ) {
					Vmatch = false;
					continue;
				}

				// If we've reached this point: V sequence matches !!
				accuBinCount += binCount;
				byteCount = ( c_n / 4 );
				binCount = (( c_n % 4 ) * 2 );

				// Compare full bytes of n with InVivo sequence
				
				if ( c_n != 0 ) {
					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = encNSeq[ m ];
						tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
						if ( tempChar0 != tempChar1 ) {
							seqMatch = false;
							break;
						}
						sh_index++;
					}

					// Compare the overflow bits of V with Invivo sequence
					tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
					tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar2 ) {
						seqMatch = false;
					}
				}

				if ( seqMatch == true ) {
					// If we've reached this point: the sequence matches
					arraylocal[2*128 + threadIdx.x] += c_DB_Full_Chew_Occur;
//					printf("Vn: i: %d, nval: %d, pairBase: %d, Vindx: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, valueAdd: %d \n", i, c_n, (pairBase / 16), Vindx, length, seqLen, g_tid, sum, c_DB_Full_Chew_Occur);
				}
			}

			if ( Vmatch == false ) {
				continue;
			}


			/////////////////////////////////////////////////////////
			//Compare InVivo Sequence to VnJ comb with D chewed
			/////////////////////////////////////////////////////////

			for ( int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++ ) {
				
				arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharV[ Vindx ] + const_d_numUniqueCharJ[ Jindx ] + c_n;
				if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] ) {
//					printf(" VnJ Length matches: i: %d, Vindx: %d, Jindx: %d, c_n: %d, Vlen: %d, Jlen: %d, length: %d, seqLen: %d\n", i, Vindx, Jindx, c_n, const_d_numUniqueCharV[ Vindx ], const_d_numUniqueCharJ[ Jindx ], length, seqLen );
					seqMatch = true;
					sh_index = 2;
					accuBinCount = 0;
					byteCount = (( const_d_numUniqueCharV[ Vindx ] ) / 4 );		// Calculates the full bytes of V
					binCount = ( const_d_numUniqueCharV[ Vindx ] % 4 ) * 2;		// Calculates the overflow bits of V
					k = const_d_V_base[Vindx];									// Starting address of V sequence

					// Compare the full bytes of V with InVivo sequence
					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = const_d_V[ k ];
						tempChar1 = iterSeq_sm[ sh_index ];
						if ( tempChar0 != tempChar1 ) {
							Vmatch = false;
							break;
						}
						sh_index++;
						k++;
					}
					
					// V does not match => get next V sequence => break J loop 
					if ( Vmatch == false ) {
						break;
					}
					
					// Compare the overflow bits of V with Invivo sequence
					tempChar0 = ((( const_d_V[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ((( iterSeq_sm[ sh_index ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar1 ) {
						Vmatch = false;
						break;
					}

					// If we've reached this point: V sequence matches !!
					accuBinCount += binCount;
					byteCount = ( c_n / 4 );
					binCount = (( c_n % 4 ) * 2 );

					if ( c_n != 0 ) {
						// Compare the full bytes of n with InVivo sequence
						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = encNSeq[ m ];
							tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));	
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							sh_index++;
						}
						
						if ( seqMatch == false ) {
							break;
						}

						// Compare the overflow bits of n with Invivo sequence
						tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							break;
						}
					}

					accuBinCount += binCount;
					if ( accuBinCount >= 8 ) {
						sh_index++;
					}
					accuBinCount %= 8;

					byteCount = (( const_d_numUniqueCharJ[ Jindx ] ) / 4 );			// Calculates the full bytes of J
					binCount = ( const_d_numUniqueCharJ[ Jindx ] % 4 ) * 2;			// Calculates the overflow bits of J
					k = const_d_J_base[ Jindx ];									// Starting address of V sequence

					// Compare the full bytes of J with InVivo sequence
					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = const_d_J[ k ];
						tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
						if ( tempChar0 != tempChar1 ) {
							seqMatch = false;
							break;
						}
						sh_index++;
						k++;
					}
					
					if ( seqMatch == false ) {
						continue;
					}

					// Compare the overflow bits of J with Invivo sequence
					tempChar0 = ((( const_d_J[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
					tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar2 ) {
						seqMatch = false;
						continue;
					}
					
					if ( seqMatch == true ) {
						arraylocal[2*128 + threadIdx.x] += c_DB_Full_Chew_Occur; 		// The sequence matches
//						printf("VnJ: i: %d, n: %d, pairBase: %d, Vindx: %d, Jindx: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, c_n, (pairBase / 16), Vindx, Jindx, length, seqLen, g_tid, sum, c_DB_Full_Chew_Occur);
					}
				}
			}

			if ( Vmatch == false ) {
				continue;
			}

			/////////////////////////////////////////////////////////
			//Compare InVivo Sequence to VnDn comb with J chewed
			/////////////////////////////////////////////////////////

			for ( int Dindx = 0; Dindx < const_numDB1; Dindx++ ) {
				
				seqMatch = true;
				arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharV[ Vindx ] + const_d_numUniqueCharDB1[ Dindx ] + c_n;
				if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] ) {
//					printf(" VnJ Length matches: i: %d, Vindx: %d, Dindx: %d, c_n: %d, Vlen: %d, Dlen: %d, length: %d, seqLen: %d\n", i, Vindx, Dindx, c_n, const_d_numUniqueCharV[ Vindx ], const_d_numUniqueCharDB1[ Dindx ], length, seqLen );
					seqMatch = true;
					sh_index = 2;
					accuBinCount = 0;
					byteCount = (( const_d_numUniqueCharV[ Vindx ] ) / 4 );		// Calculates the full bytes of V
					binCount = ( const_d_numUniqueCharV[ Vindx ] % 4 ) * 2;		// Calculates the overflow bits of V
					k = const_d_V_base[Vindx];									// Starting address of V sequence

					// Compare the full bytes of V with InVivo sequence
					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = const_d_V[ k ];
						tempChar1 = iterSeq_sm[ sh_index ];
						if ( tempChar0 != tempChar1 ) {
							Vmatch = false;
							break;
						}
						sh_index++;
						k++;
					}
					
					// V does not match => get next V sequence => break J loop 
					if ( Vmatch == false ) {
						break;
					}
					
					// Compare the overflow bits of V with Invivo sequence
					tempChar0 = ((( const_d_V[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ((( iterSeq_sm[ sh_index ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar1 ) {
						Vmatch = false;
						break;
					}

//					if ( i == 157 ) {
//						printf("VnD Byte match: i: %d, pairBase: %d, c_n: %d, Vindx: %d, Dindx: %d, Vlen: %d, Dlen: %d, length: %d, seqLen: %d, g_tid: %d\n", i, (pairBase / 16), c_n, Vindx, Dindx, const_d_numUniqueCharV[ Vindx ], const_d_numUniqueCharDB1[ Dindx ], length, seqLen, g_tid);
//					}

					accuBinCount += binCount;
					
					for ( int nlen = 0; nlen < ( c_n + 1 ); nlen++ ) {

						seqMatch = true;
						tempAccuBinCount = accuBinCount;
						tempShIndex = sh_index;
						byteCount = ( nlen / 4 );
						binCount = (( nlen % 4 ) * 2 );

						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = encNSeq[ m ];
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));	
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
						}

						if ( seqMatch == false ) {
							continue;
						}
						
						// Compare the overflow bits of n with Invivo sequence
						tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							continue;
						}
						
						tempAccuBinCount += binCount;
						if ( tempAccuBinCount >= 8 ) {
							tempShIndex++;
						}
						tempAccuBinCount %= 8;

						byteCount = ( const_d_numUniqueCharDB1[ Dindx ] / 4 );
						binCount = (( const_d_numUniqueCharDB1[ Dindx ] % 4 ) * 2 );
						k = const_d_DB1_base[ Dindx ];									// Starting address of V sequence

						// Compare the full bytes of D with InVivo sequence
						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = const_d_DB1[ k ];
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
							k++;
						}

						if ( seqMatch == false ) {
							continue;
						}
						
						// Compare the overflow bits of D with Invivo sequence
						tempChar0 = ((( const_d_DB1[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							continue;
						}

						tempAccuBinCount += binCount;
						if ( tempAccuBinCount >= 8 ) {
							tempShIndex++;
						}
						tempAccuBinCount %= 8;

//						printf("VnD Byte match: i: %d, nlen: %d, c_n: %d, Vindx: %d, Dindx: %d, g_tid: %d\n", i, nlen, c_n, Vindx, Dindx, g_tid );
//						Worked till here!

						byteCount = ( nlen / 4 );
						binCount = (( nlen % 4 ) * 2 );

						for ( int m = 0; m < (( c_n - nlen ) / 4 ); m++ ) {
							tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
							byteCount++;
						}

						if ( seqMatch == false ) {
							continue;
						}

						tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
//						if ( i == 157 ) {
//							printf("i: %d, nlen: %d, c_n: %d, Vindx: %d, Dindx: %d, g_tid: %d \n\t byteCount: %d, binCount: %d, tempAccuBinCount: %d, tempShIndex: %d \n\t encNSeq: %d, tempChar0: %d \n\t iterSeq_sm: %d %d, tempChar1: %d, tempChar2: %d \n", i, nlen, c_n, Vindx, Dindx, g_tid, byteCount, binCount, tempAccuBinCount, tempShIndex, encNSeq[ byteCount ], tempChar0, iterSeq_sm[ tempShIndex ], iterSeq_sm[ tempShIndex + 1 ], tempChar1, tempChar2);
//						}
						binCount = ((( c_n - nlen ) % 4 ) * 2 );
						tempChar0 = ((( tempChar0 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar1 ) {
							seqMatch = false;
							continue;
						}

						if ( seqMatch == true )	{
							arraylocal[2*128 + threadIdx.x] += const_d_numOccurrenceDB1[ Dindx ];
//							printf("VnDn: i: %d,pairBase: %d, Vindx: %d, Dindx: %d, nlen: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, (pairBase / 16), Vindx, Dindx, nlen, length, seqLen, g_tid, sum, const_d_numOccurrenceDB1[Dindx]);
						}
					}

				}

				/////////////////////////////////////////////////////////
				//Compare InVivo Sequence to VnDnJ comb with no chewback
				/////////////////////////////////////////////////////////

				for ( int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++ ) {	

					seqMatch = true;
					arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharV[Vindx] + const_d_numUniqueCharJ[Jindx] + const_d_numUniqueCharDB1[Dindx] + c_n;
					if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] )
					{	
						sh_index = 2;													// InVivo sequence starts from iterSeq_sm[2]
						accuBinCount = 0;
						byteCount = (( const_d_numUniqueCharV[ Vindx ] ) / 4 );			// Calculates the full bytes of V
						binCount = ( const_d_numUniqueCharV[ Vindx ] % 4 ) * 2;			// Calculates the overflow bits of V
						k = const_d_V_base[Vindx];										// Starting address of V sequence

						// Compare the full bytes of V with InVivo sequence
						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = const_d_V[ k ];
							tempChar1 = iterSeq_sm[ sh_index ];
							if ( tempChar0 != tempChar1 ) {
								Vmatch = false;
								break;
							}
							sh_index++;
							k++;
						}
						
						// V does not match => get next V sequence => break J loop 
						if ( Vmatch == false ) {
							break;
						}
						
						// Compare the overflow bits of V with Invivo sequence
						tempChar0 = ((( const_d_V[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ((( iterSeq_sm[ sh_index ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar1 ) {
							Vmatch = false;
							break;
						}
						
//						printf("VnDnJ_Vmatch: i: %d, nval: %d, pairBase: %d, Vindx: %d, Dindx: %d, Jindx: %d, seqLen: %d, Vlen: %d, Dlen: %d, Jlen: %d, g_tid: %d \n", i, c_n, (pairBase / 16), Vindx, Dindx, Jindx, seqLen, const_d_numUniqueCharV[Vindx], const_d_numUniqueCharDB1[Dindx], const_d_numUniqueCharJ[Jindx], g_tid);

						accuBinCount += binCount;
						
						for ( int nlen = 0; nlen < ( c_n + 1 ); nlen++ ) {

							seqMatch = true;
							tempAccuBinCount = accuBinCount;
							tempShIndex = sh_index;
							byteCount = ( nlen / 4 );
							binCount = (( nlen % 4 ) * 2 );

							for ( int m = 0; m < byteCount; m++ ) {
								tempChar0 = encNSeq[ m ];
								tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));	
								if ( tempChar1 != tempChar0 ) {
									seqMatch = false;
									break;
								}
								tempShIndex++;
							}	

							if ( seqMatch == false ) {
								continue;
							}
						
							// Compare the overflow bits of n with Invivo sequence
							tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							if ( tempChar0 != tempChar2 ) {
								seqMatch = false;
								continue;
							}
						
							tempAccuBinCount += binCount;
							if ( tempAccuBinCount >= 8 ) {
								tempShIndex++;
							}
							tempAccuBinCount %= 8;

							byteCount = ( const_d_numUniqueCharDB1[ Dindx ] / 4 );
							binCount = (( const_d_numUniqueCharDB1[ Dindx ] % 4 ) * 2 );
							k = const_d_DB1_base[ Dindx ];									// Starting address of V sequences

							// Compare the full bytes of D with InVivo sequence
							for ( int m = 0; m < byteCount; m++ ) {
								tempChar0 = const_d_DB1[ k ];
								tempChar2 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
								if ( tempChar0 != tempChar2 ) {
									seqMatch = false;
									break;
								}
								tempShIndex++;
								k++;
							}

							if ( seqMatch == false ) {
								continue;
							}
						
							// Compare the overflow bits of D with Invivo sequence
							tempChar0 = ((( const_d_DB1[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							if ( tempChar0 != tempChar2 ) {
								seqMatch = false;
								continue;
							}

							tempAccuBinCount += binCount;
							if ( tempAccuBinCount >= 8 ) {
								tempShIndex++;
							}
							tempAccuBinCount %= 8;

							byteCount = ( nlen / 4 );
							binCount = (( nlen % 4 ) * 2 );	

							for ( int m = 0; m < (( c_n - nlen ) / 4 ); m++ ) {
								tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
								tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
								if ( tempChar0 != tempChar1 ) {
									seqMatch = false;
									break;
								}
								tempShIndex++;
								byteCount++;
							}

							if ( seqMatch == false ) {
								continue;
							}

							tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							
							binCount = ((( c_n - nlen ) % 4 ) * 2 );
							tempChar0 = ((( tempChar0 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								continue;
							}

							tempAccuBinCount += binCount;
							if ( tempAccuBinCount >= 8 ) {
								tempShIndex++;
							}
							tempAccuBinCount %= 8;

							byteCount = (( const_d_numUniqueCharJ[ Jindx ] ) / 4 );			// Calculates the full bytes of J
							binCount = ( const_d_numUniqueCharJ[ Jindx ] % 4 ) * 2;			// Calculates the overflow bits of J
							k = const_d_J_base[ Jindx ];									// Starting address of V sequence

							for ( int m = 0; m < byteCount; m++ ) {
								tempChar0 = const_d_J[ k ];
								tempChar2 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
								if ( tempChar0 != tempChar2 ) {
									seqMatch = false;
									break;
								}
								tempShIndex++;
								k++;
							}

							if ( seqMatch == false ) {
								break;
							}

							tempChar0 = ((( const_d_J[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							if ( tempChar0 != tempChar2 ) {
								seqMatch = false;
								break;
							}

							if ( seqMatch == true ) {
								arraylocal[2*128 + threadIdx.x] += const_d_numOccurrenceDB1[ Dindx ];
//								printf("VnDnJ: i: %d, n: %d, pairBase: %d, Vindx: %d, Jindx: %d, Dindx: %d, nlen: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, c_n, (pairBase / 16), Vindx, Jindx, Dindx, nlen, length, seqLen, g_tid, sum, const_d_numOccurrenceDB1[Dindx]); 
							}
						}
					}
				}

				if ( Vmatch == false ) {
					break;
				}

			}

			if ( Vmatch == false ) {
				continue;
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		//Compare our InVivo Sequences containing 'nDn' = nDn, nDnJ
		//////////////////////////////////////////////////////////////////////////////////

		// Get all the D sequences:
		for ( int Dindx = 0; Dindx < const_numDB1; Dindx++ ) {

			arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharDB1[Dindx] + c_n;
			if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] )
			{	

				//int tempAccuBinCount = 0;
				//int tempShIndex = 0;

				for ( int nlen = 0; nlen < ( c_n + 1 ); nlen++ ) {

					seqMatch = true;
					tempAccuBinCount = 0;
					tempShIndex = 2;
					byteCount = ( nlen / 4 );
					binCount = (( nlen % 4 ) * 2 );

					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = encNSeq[ m ];
						tempChar2 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));	
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							break;
						}
						tempShIndex++;
					}	

					if ( seqMatch == false ) {
						continue;
					}
						
					// Compare the overflow bits of n with Invivo sequence
					tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
					tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar2 ) {
						seqMatch = false;
						continue;
					}
						
					tempAccuBinCount += binCount;
					if ( tempAccuBinCount >= 8 ) {
						tempShIndex++;
					}
					tempAccuBinCount %= 8;

					byteCount = ( const_d_numUniqueCharDB1[ Dindx ] / 4 );
					binCount = (( const_d_numUniqueCharDB1[ Dindx ] % 4 ) * 2 );
					k = const_d_DB1_base[ Dindx ];									// Starting address of D sequences

					// Compare the full bytes of D with InVivo sequence
					for ( int m = 0; m < byteCount; m++ ) {
						tempChar0 = const_d_DB1[ k ];
						tempChar2 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							break;
						}
						tempShIndex++;
						k++;
					}

					if ( seqMatch == false ) {
						continue;
					}
						
					// Compare the overflow bits of D with Invivo sequence
					tempChar0 = ((( const_d_DB1[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
					tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar2 ) {
						seqMatch = false;
						continue;
					}

					tempAccuBinCount += binCount;
					if ( tempAccuBinCount >= 8 ) {
						tempShIndex++;
					}
					tempAccuBinCount %= 8;

					byteCount = ( nlen / 4 );
					binCount = (( nlen % 4 ) * 2 );	

					for ( int m = 0; m < (( c_n - nlen ) / 4 ); m++ ) {
						tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						if ( tempChar0 != tempChar1 ) {
							seqMatch = false;
							break;
						}
						tempShIndex++;
						byteCount++;
					}

					if ( seqMatch == false ) {
						continue;
					}

					tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							
					binCount = ((( c_n - nlen ) % 4 ) * 2 );
					tempChar0 = ((( tempChar0 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					tempChar1 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
					if ( tempChar0 != tempChar1 ) {
						seqMatch = false;
						continue;
					}

					if ( seqMatch == true ) {
						arraylocal[2*128 + threadIdx.x] += const_d_numOccurrenceDB1[ Dindx ];
//						printf("nDn: i: %d, n: %d, pairBase: %d, Dindx: %d, nlen: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, c_n, (pairBase / 16), Dindx, nlen, length, seqLen, g_tid, sum, const_d_numOccurrenceDB1[ Dindx ]);
					}
				}
			}

			/////////////////////////////////////////////////////////
			//Compare InVivo Sequence to nDnJ comb with V chewed
			/////////////////////////////////////////////////////////

			for ( int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++ ) {
//				printf("nDnJ: i: %d, nval: %d, pairBase: %d,Dindx: %d, Jindx: %d, seqLen: %d, Dlen: %d, Jlen: %d, g_tid: %d \n", i, c_n, (pairBase / 16), Dindx, Jindx, seqLen, const_d_numUniqueCharDB1[Dindx], const_d_numUniqueCharJ[Jindx], g_tid);
				arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharJ[Jindx] + const_d_numUniqueCharDB1[Dindx] + c_n;
				if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] )
				{	

					//int tempAccuBinCount = 0;
					//int tempShIndex = 0;

					for ( int nlen = 0; nlen < ( c_n + 1 ); nlen++ ) {

						seqMatch = true;
						tempAccuBinCount = 0;
						tempShIndex = 2;
						byteCount = ( nlen / 4 );
						binCount = (( nlen % 4 ) * 2 );

						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = encNSeq[ m ];
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));	
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
						}	

						if ( seqMatch == false ) {
							continue;
						}
					
						// Compare the overflow bits of n with Invivo sequence
						tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							continue;
						}
						
						tempAccuBinCount += binCount;
						if ( tempAccuBinCount >= 8 ) {
							tempShIndex++;
						}
						tempAccuBinCount %= 8;

						byteCount = ( const_d_numUniqueCharDB1[ Dindx ] / 4 );
						binCount = (( const_d_numUniqueCharDB1[ Dindx ] % 4 ) * 2 );
						k = const_d_DB1_base[ Dindx ];									// Starting address of V sequences

						// Compare the full bytes of D with InVivo sequence
						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = const_d_DB1[ k ];
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
							k++;
						}

						if ( seqMatch == false ) {
							continue;
						}
						
						// Compare the overflow bits of D with Invivo sequence
						tempChar0 = ((( const_d_DB1[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							continue;
						}

						tempAccuBinCount += binCount;
						if ( tempAccuBinCount >= 8 ) {
							tempShIndex++;
						}
						tempAccuBinCount %= 8;

						byteCount = ( nlen / 4 );
						binCount = (( nlen % 4 ) * 2 );	

						for ( int m = 0; m < (( c_n - nlen ) / 4 ); m++ ) {
							tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
							tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							if ( tempChar0 != tempChar1 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
							byteCount++;
						}

						if ( seqMatch == false ) {
							continue;
						}

						tempChar0 = ( encNSeq[ byteCount ] << binCount ) | ((( encNSeq[ byteCount + 1 ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
								
						binCount = ((( c_n - nlen ) % 4 ) * 2 );
						tempChar0 = ((( tempChar0 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						tempChar1 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar1 ) {
							seqMatch = false;
							continue;
						}

						tempAccuBinCount += binCount;
						if ( tempAccuBinCount >= 8 ) {
							tempShIndex++;
						}
						tempAccuBinCount %= 8;

						byteCount = (( const_d_numUniqueCharJ[ Jindx ] ) / 4 );			// Calculates the full bytes of J
						binCount = ( const_d_numUniqueCharJ[ Jindx ] % 4 ) * 2;		// Calculates the overflow bits of J
						k = const_d_J_base[ Jindx ];									// Starting address of V sequence

						for ( int m = 0; m < byteCount; m++ ) {
							tempChar0 = const_d_J[ k ];
							tempChar2 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | ((( iterSeq_sm[ tempShIndex + 1 ] >> 2 ) & 0x3F ) >> ( 6 - tempAccuBinCount ));
							if ( tempChar0 != tempChar2 ) {
								seqMatch = false;
								break;
							}
							tempShIndex++;
							k++;
						}

						if ( seqMatch == false ) {
							break;
						}

						tempChar0 = ((( const_d_J[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount )) ;
						tempChar1 = ( iterSeq_sm[ tempShIndex ] << tempAccuBinCount ) | (( iterSeq_sm[ tempShIndex + 1] >> 2 ) & 0x3F >> ( 6 - tempAccuBinCount )); 
						tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
						if ( tempChar0 != tempChar2 ) {
							seqMatch = false;
							break;
						}

						if ( seqMatch == true ) {
							arraylocal[2*128 + threadIdx.x] += const_d_numOccurrenceDB1[ Dindx ];
//							printf("nDnJ: i: %d, n: %d, pairBase: %d, Jindx: %d, Dindx: %d, nlen: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, c_n, (pairBase / 16), Jindx, Dindx, nlen, length, seqLen, g_tid, sum, const_d_numOccurrenceDB1[ Dindx ]);
						}
					}
				}
			}
		}

		/////////////////////////////////////////////////////////
		//Compare InVivo Sequence to n comb with V, D and J chewed
		/////////////////////////////////////////////////////////

//		printf("n: i: %d, nval: %d, pairBase: %d, seqLen: %d, g_tid: %d \n", i, c_n, (pairBase / 16), seqLen, g_tid);
		if ( c_n == arraylocal[6*128 + threadIdx.x] )
		{	
			sh_index = 2;
			accuBinCount = 0;

			byteCount = ( c_n / 4 );
			binCount = (( c_n % 4 ) * 2 );

			// Compare full bytes of n with InVivo sequence
			
			for ( int m = 0; m < byteCount; m++ ) {
				tempChar0 = encNSeq[ m ];
				tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));	
				if ( tempChar1 != tempChar0 ) {
					seqMatch = false;
					break;
				}
				sh_index++;
			}

			if ( seqMatch == false ) {
				break;
			}

			// Compare the overflow bits of V with Invivo sequence
			tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
			tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
			tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
			if ( tempChar0 != tempChar2 ) {
				seqMatch = false;
			}

			if ( seqMatch == true ) {
				arraylocal[2*128 + threadIdx.x] += c_DB_Full_Chew_Occur;
//				printf("n: i: %d, pairBase: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, (pairBase / 16), length, seqLen, g_tid, sum, c_DB_Full_Chew_Occur);
			}
		}

		/////////////////////////////////////////////////////////
		//Compare InVivo Sequence to nJ comb with V and D chewed
		/////////////////////////////////////////////////////////
		for ( int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++ ) {

//			printf("nJ: i: %d, nval: %d, pairBase: %d, Jindx: %d, seqLen: %d, Jlen: %d, g_tid: %d \n", i, c_n, (pairBase / 16), Jindx, seqLen, const_d_numUniqueCharJ[Jindx], g_tid);
			arraylocal[8*128 + threadIdx.x] = const_d_numUniqueCharJ[Jindx] + c_n;
			if ( arraylocal[8*128 + threadIdx.x] == arraylocal[6*128 + threadIdx.x] )
			{	
				sh_index = 2;
				accuBinCount = 0;

				byteCount = ( c_n / 4 );
				binCount = (( c_n % 4 ) * 2 );

				// Compare full bytes of n with InVivo sequence
				for ( int m = 0; m < byteCount; m++ ) {
					tempChar0 = encNSeq[ m ];
					tempChar2 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));	
					if ( tempChar2 != tempChar0 ) {
						seqMatch = false;
						break;
					}
					sh_index++;
				}

				if ( seqMatch == false ) {
					break;
				}

				// Compare the overflow bits of V with Invivo sequence
				tempChar0 = ((( encNSeq[ byteCount ] >> 2 ) & 0x3F ) >> ( 6 - binCount ));
				tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
				tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
				if ( tempChar0 != tempChar2 ) {
					seqMatch = false;
					break;
				}

				accuBinCount += binCount;
				if ( accuBinCount >= 8 ) {
					sh_index++;
				}
				accuBinCount %= 8;

				byteCount = (( const_d_numUniqueCharJ[ Jindx ] ) / 4 );			// Calculates the full bytes of J
				binCount = ( const_d_numUniqueCharJ[ Jindx ] % 4 ) * 2;			// Calculates the overflow bits of J
				k = const_d_J_base[ Jindx ];									// Starting address of V sequence

				// Compare the full bytes of J with InVivo sequence
				for ( int m = 0; m < byteCount; m++ ) {
					tempChar0 = const_d_J[ k ];
					tempChar2 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | ((( iterSeq_sm[ sh_index + 1 ] >> 2 ) & 0x3F ) >> ( 6 - accuBinCount ));
					if ( tempChar0 != tempChar2 ) {
						seqMatch = false;
						break;
					}
					sh_index++;
					k++;
				}

				if ( seqMatch == false ) {
					continue;
				}

				tempChar0 = ((( const_d_J[ k ] >> 2 ) & 0x3F ) >> ( 6 - binCount )) ;
				tempChar1 = ( iterSeq_sm[ sh_index ] << accuBinCount ) | (( iterSeq_sm[ sh_index + 1] >> 2 ) & 0x3F >> ( 6 - accuBinCount )); 
				tempChar2 = ((( tempChar1 >> 2 ) & 0x3F ) >> ( 6 - binCount ));
				if ( tempChar0 != tempChar2 ) {
					seqMatch = false;
					continue;
				}

				if ( seqMatch == true ) {
					arraylocal[2*128 + threadIdx.x] += c_DB_Full_Chew_Occur;
//					printf("nJ: i: %d, Jnum: %d, n: %d, pairBase: %d, Jindx: %d, length: %d, seqLen: %d, g_tid: %d, sum: %d, addVal: %d \n", i, Jnum, c_n, (pairBase / 16), Jindx, length, seqLen, g_tid, sum, c_DB_Full_Chew_Occur);
				}
			}
		}

		//-------------------------------------------------------------------------------------------------------
		//If only 1 thread-block, then we can write results to RAM using just InVivo sequence number
		//-------------------------------------------------------------------------------------------------------
		if( blockDim.x == 1 ) {		//if there is only 1 thread per block, just use i as global memory index. No need for reduction
			d_Results[ i ] = arraylocal[2*128 + threadIdx.x];
		}


		//-------------------------------------------------------------------------------------------------------
		//					Perform Reduction of Results if more than 1 thread-block				  
		//-------------------------------------------------------------------------------------------------------
		//reduction for current InVivo sequence in shared memory
		if( blockDim.x > 1 ) {
	
			result_sm[ threadIdx.x ] = arraylocal[2*128 + threadIdx.x];			//write a threads sum to the shared memory
			__syncthreads();						//make sure all sums have been written before proceeding

			int half = blockDim.x / 2;

			while( 1 ) {								//how many reductions we need
				if( threadIdx.x < half ) {				//only certain threads perform reduction
					result_sm[ threadIdx.x ] += result_sm[ threadIdx.x + half ];
				}
				__syncthreads();
				if( half == 1 ) break;
				half = half / 2;
			}

			__syncthreads();

			//write results to the global memory. Each thread-block writes 1 result for each InVivo Sequence i
			if( threadIdx.x == 0 ){					//we need only 1 thread in the thread block to write its result
				d_Results[ i * gridDim.x + blockIdx.x ] = result_sm[ 0 ];		//write our consolidated result into the global memory
			}

		}	//end result reduction
	}		//end iterating through InVivo Sequences
	
	return;
} 			//kernel done

#endif // #ifndef _TNT_KERNEL_H_


//	printf("Block Dim: %d, Block Idx: %d, Thread Idx: %d, g_tid: %d, encNSeq: %d,%d,%d\n", blockDim.x, blockIdx.x, threadIdx.x, g_tid, encNSeq[0], encNSeq[1], encNSeq[2]);

//	if ( g_tid == 0 ) {
//		printf("num_Seqs: %d, pairBase: %d\n", num_Seqs, pairBase);
//	}

//	if ( threadIdx.x == 1 ) {
//		printf("i: %d, threadIdx: %d, iterseq_sm: %d, gl_index: %d, d_InVivo_cp64: %d\n", i, threadIdx.x, iterSeq_sm[ threadIdx.x ], gl_index, d_InVivo_cp64[ gl_index ]);
//	}

//	if ( g_tid == 12 ) {
//		printf("i: %d, iterSeq_sm: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", i, iterSeq_sm[0], iterSeq_sm[1], iterSeq_sm[2], iterSeq_sm[3], iterSeq_sm[4], iterSeq_sm[5], iterSeq_sm[6], iterSeq_sm[7], iterSeq_sm[8], iterSeq_sm[9], iterSeq_sm[10], iterSeq_sm[11], iterSeq_sm[12], iterSeq_sm[13], iterSeq_sm[14], iterSeq_sm[15]);
//	}

//	__syncthreads();

//	if ( g_tid == 0 ) {
//		printf("i: %d, getChar: %d, seqLen: %d \n", i, getChar, seqLen);
//	}

//	if ( g_tid == 4 ) {
//		printf("getChar: %d, const_d_V[k]: %d \n", getChar, const_d_V[k]);
//	}

//	if ( g_tid == 4 ) {
//		printf("Byte V matches: i: %d, Vnum: %d, pairBase: %d, Vindx: %d, length: %d, seqLen: %d \n", i, Vnum, (pairBase / 16), Vindx, length, seqLen);
//	}

//	if ( g_tid == 4 ) {
//		printf("const_d_V: %d, iterSeq_sm: %d, V_tempChar0: %d, V_tempChar1: %d \n", const_d_V[k], iterSeq_sm[sh_index], tempChar0, tempChar1);
//	}

//	if ( g_tid == 4 ) {
//		printf("binCount V matches: i: %d, Vnum: %d, pairBase: %d, Vindx: %d, length: %d, seqLen: %d \n", i, Vnum, (pairBase / 16), Vindx, length, seqLen);
//	}

//	if ( seqMatch == true ) {
//		printf("Byte n matches: i: %d, Vnum: %d, pairBase: %d, Vindx: %d, length: %d, seqLen: %d, g_tid: %d \n", i, Vnum, (pairBase / 16), Vindx, length, seqLen, g_tid);
//	}

//	if ( seqMatch == true ) {
//		printf("n_tempChar0: %d, n_tempChar2: %d, g_tid: %d \n", tempChar0, tempChar2, g_tid);
//		printf("binCount n matches: i: %d, Vnum: %d, pairBase: %d, Vindx: %d, length: %d, seqLen: %d, g_tid: %d \n", i, Vnum, (pairBase / 16), Vindx, length, seqLen, g_tid);
//	}

