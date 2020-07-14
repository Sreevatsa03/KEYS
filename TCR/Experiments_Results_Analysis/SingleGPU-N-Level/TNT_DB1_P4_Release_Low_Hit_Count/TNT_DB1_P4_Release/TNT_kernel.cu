#ifndef _TNT_KERNEL_H_
#define _TNT_KERNEL_H_

#include <stdio.h>

//DB information
__constant__ int const_numDB1 = 169;				//total number of DB1 chewbacks
__constant__ char const_d_DB1[1448];				//constant memory allocation for DB chewbacks minus full chew back
__constant__ int const_d_DB1_base[169];				//constant memory contains location of each starting sequence in d_DB1 
__constant__ int const_d_numOccurrenceDB1[169];		//Number of ways for particular DB1 chewback 
__constant__ int const_d_numUniqueCharDB1[169];		//number of characters in a unique occurence of DB1
//V information
__constant__ int const_numV = 342;					//total number of V sequences in all V files
__constant__ char const_d_V[3107];					//holds all V chewback sequences
__constant__ int const_d_V_base[342];				//contains the starting index of each V sequence
__constant__ int const_d_numUniqueCharV[342];		//number of characters in a unique occurence of V
//J information
__constant__ int const_numJ = 271;					//total number of J sequences in all J files
__constant__ char const_d_J[3210];					//holds all J chewback sequences
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


//char d_chewArrV[1500][30];
//char d_chewArrJ[1500][30];
//int d_counterIndex = 0;
__device__ int d_counterIndex;




/////////////////////////////////////////////////
//kernel for 64 threads or less
/////////////////////////////////////////////////
__global__ void
TNT_kernel_InVivo64(unsigned int* d_Results, char* d_InVivo_cp64, char* d_chewArrV, char* d_chewArrJ)	
{

	volatile __shared__ char iterSeq_sm[64]; //the thread block size we will use for this kernel is 64
	volatile __shared__ int result_sm[512];  //the max thread-block size

	//The four possible bases
	char base[4] = {'A', 'T', 'G', 'C'};

	char nSeq[12];					//will hold a single n combination
	int Vnum = c_Vnum;				//current V file
	int Jnum = c_Jnum;				//current J file
	int sh_index;					//used as a shared memory index
	int sum;						//holds an iterative sum for result
	char tmpChar;					//used to temporarily hold a character


	//obtain a unique global index for each thread in the grid
	unsigned int g_tid = blockIdx.x*blockDim.x + threadIdx.x;

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

	//get the number of InVivo VJ sequences we need to go through
	int num_Seqs = const_d_VJ_Pairs[Vnum*12 + Jnum];				//multiply by 12. Number of J files.

	//int whichSeq; //which sequence is our current thread-block working on in the scope of current VJ
	int seqLen;	  //length of our current sequence
	int pairBase = const_VJ_Pair_Base[Vnum*12 + Jnum] * 64; //The base address for a given VJ pair


	//iterate through all InVivo combinations for current VJ pair
	for(int i = 0; i < num_Seqs; i++){

		result_sm[threadIdx.x] = 0;
		sum = 0; //reset our result

		__syncthreads();

		//store an InVivo combination into the shared memory "iterResults_sm[]"
		if(blockDim.x < 64){								//iter through VJ seq if block dim < 64. There's only 1 block
			for(int j = 0; j < (64 / blockDim.x); j++){		//iterations = sequence allocation / block size
				int k = j*blockDim.x+threadIdx.x;			//create local SM index
				int gl_index = (pairBase + i*64) + k;		//create global memory index
				iterSeq_sm[k] = d_InVivo_cp64[gl_index];    //read the current InVivo sequence from the global memory
			}
		}
		else{	//only threads < 64 will read inVivo data
			if(threadIdx.x < 64){
				int gl_index = (pairBase + i*64) + threadIdx.x;	  //create global memory index
				iterSeq_sm[threadIdx.x] = d_InVivo_cp64[gl_index];//read the current InVivo sequence from the global memory
			}
		}



		//if(blockDim.x > 1) 
			__syncthreads();

		//get the length of current sequence for all threads in current thread-block
		seqLen = (int)iterSeq_sm[2];

		//set our shared memory index to the base of the sequence characters in shared array
		sh_index = 3;

		int n_p1 = c_n + 1;
		int n_cnt, k;
		int length;			//length of each sequence we generate
		bool Vmatch = 1;	//Is there a V sequence match?
		bool seqMatch = 1;	//Is the entire sequence a match?


		//////////////////////////////////////////////////////////////////////////////////
		//First compare our InVivo Sequences to VnDnJ combinations with no full chewbacks
		//////////////////////////////////////////////////////////////////////////////////
		for(int Vindx = c_V_Begin; Vindx < c_V_End; Vindx++){			//go through relevent V sequences
		
			Vmatch = 1;	 //assume V is a match before we check it
			seqMatch = 1;

			/////////////////////////////////////////////////////////
			//Compare InVivo Sequence to Vn comb with D and J chewed 
			/////////////////////////////////////////////////////////
			length = const_d_numUniqueCharV[Vindx] + c_n;
			n_cnt = c_n;
					
			//check to see if sequence to create is the same length as the InVivo sequence
			if(seqLen == length){
				sh_index = 3;											//reset our shared memory index

				k = const_d_V_base[Vindx];									//starting address of V sequence
				for(int m = 0; m < const_d_numUniqueCharV[Vindx]; m++){		//go through each character in current V sequence
					tmpChar = const_d_V[k];									//load a V character into a temp variable
					if(tmpChar != iterSeq_sm[sh_index]){Vmatch = 0; break;} //End V comparisons
					sh_index++;												//increment shared memory index
					k++;													//increment for next character
				}
				
				if(Vmatch == 0) continue; //check the next V sequence

				if(c_n != 0){
					//add n combination
					for(int m = 0; m < n_cnt; m++){
						tmpChar = nSeq[m]; 
						if(tmpChar != iterSeq_sm[sh_index]){seqMatch = 0; break;}	//jump to next iteration if a character does not match
						sh_index++;
					}
				}

				if (seqMatch) 
					sum += c_DB_Full_Chew_Occur; //sum++; //if we've made it this far, the sequences match.
			}

			if(Vmatch == 0) continue;


			///////////////////////////////////
			//check D and J combinations
			///////////////////////////////////
			for(int Dindx = 0; Dindx < const_numDB1 && Vmatch; Dindx++){			//go through all D sequences

				/////////////////////////////////////////
				//compare VnDn with J fully chewed back 
				/////////////////////////////////////////
				length = const_d_numUniqueCharV[Vindx] + const_d_numUniqueCharDB1[Dindx] + c_n;
				n_cnt = c_n;
					
				//check to see if sequence to create is the same length as the InVivo sequence
				if(seqLen == length){

					for(int j = 0; j < n_p1 && Vmatch; j++){								//go through each n addition (n + 1)
						
						seqMatch = 1;	//Assume initially a sequence match
						sh_index = 3;	//reset our shared memory index
	
						k = const_d_V_base[Vindx];								//starting address of V sequence
						for(int m = 0; m < const_d_numUniqueCharV[Vindx]; m++){	//go through each character in current V sequence
							tmpChar = const_d_V[k];								//load a V character into a temp variable
							if(tmpChar != iterSeq_sm[sh_index]){Vmatch = 0; break;}//End V comparisons
							sh_index++;											//increment shared memory index
							k++;												//increment for next character
						}

						if(Vmatch == 0) continue;
				
						if(c_n != 0){
							//add n combination
							for(int m = 0; m < n_cnt; m++){
								tmpChar = nSeq[m]; 
								if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
								sh_index++;
							}
						}

						if(seqMatch == 0) continue;

						//glue current D sequence
						k = const_d_DB1_base[Dindx];							//starting address of D sequence
						for(int m = 0; m < const_d_numUniqueCharDB1[Dindx]; m++){	//go through each character in current V sequence
							tmpChar = const_d_DB1[k];							//store V character in shared memory
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
							sh_index++;											//increment shared memory index
							k++;												//increment for next character
						}

						if(seqMatch == 0) continue;

						if(c_n != 0){
							//add n combination
							for(int m = n_cnt; m < c_n; m++){
								tmpChar = nSeq[m]; 
								if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
								sh_index++;
							}
						}

						if(seqMatch == 0) continue;

						n_cnt--; 
						 
						sum += const_d_numOccurrenceDB1[Dindx]; //if we've made it this far, the sequences match.

	    			} //end iterating through n sequences
				} //end checking VnDn comparisons


				///////////////////////////////////
				//check D J combinations
				///////////////////////////////////
				for(int Jindx = c_J_Begin; Jindx < c_J_End && Vmatch; Jindx++){	//go through relevent J sequences

					length = const_d_numUniqueCharV[Vindx] + const_d_numUniqueCharJ[Jindx] + const_d_numUniqueCharDB1[Dindx] + c_n;
					n_cnt = c_n;

					//check to see if sequence to create is the same length as the InVivo sequence
					if(seqLen != length) continue;
	
					////////////////////////////////////////////////////
					//begin generating sequences with no full chewbacks
					////////////////////////////////////////////////////
					for(int j = 0; j < n_p1; j++){								//go through each n addition (n + 1)

						seqMatch = 1;	//Assume initially a sequence match
						sh_index = 3;	//reset our shared memory index

						k = const_d_V_base[Vindx];												//starting address of V sequence
						for(int m = 0; m < const_d_numUniqueCharV[Vindx]; m++){					//go through each character in current V sequence
							tmpChar = const_d_V[k];												//load a V character into a temp variable
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; Vmatch = 0; break;};	//End V comparisons
								sh_index++;														//increment shared memory index
								k++;															//increment for next character
						}

						if(Vmatch == 0)break; //exit current v comparisons if V is not a match
						
						if(c_n != 0){
							//add n combination
							for(int m = 0; m < n_cnt; m++){
								tmpChar = nSeq[m]; 
								if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
								sh_index++;
							}
						}
						
						if(seqMatch == 0) continue;

						//glue current D sequence
						k = const_d_DB1_base[Dindx];							//starting address of D sequence
						for(int m = 0; m < const_d_numUniqueCharDB1[Dindx]; m++){	//go through each character in current V sequence
							tmpChar = const_d_DB1[k];							//store V character in shared memory
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
							sh_index++;											//increment shared memory index
							k++;												//increment for next character
						}

						if(seqMatch == 0) continue;

						if(c_n != 0){
							//add n combination
							for(int m = n_cnt; m < c_n; m++){
								tmpChar = nSeq[m]; 
								if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	//jump to next iteration if a character does not match
								sh_index++;
							}
						}

						if(seqMatch == 0) continue;

						//glue current J sequence
						k = const_d_J_base[Jindx];								//starting address of V sequence
						for(int m = 0; m < const_d_numUniqueCharJ[Jindx]; m++){	//go through each character in current V sequence
							tmpChar = const_d_J[k];
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; seqMatch = 0; break;}	    //jump to next iteration if a character does not match
							sh_index++;											//increment shared memory index
							k++;												//increment for next character
						}

						if(seqMatch == 0) continue;

						n_cnt--; 
 
						sum += const_d_numOccurrenceDB1[Dindx]; //if we've made it this far, the sequences match.

	    			} //end iterating through n sequences
				} //end iterating through j sequences
			} //end iterating through d sequences
		} //end iterating through v sequences
	

		//---------------------------------------------------------------------------------
		//Compare our InVivo Sequences to VnJ with D full chewback
		//---------------------------------------------------------------------------------
		for(int Vindx = c_V_Begin; Vindx < c_V_End; Vindx++){			//go through relevent V sequences
			Vmatch = 1;
				for(int Jindx = c_J_Begin; Jindx < c_J_End && Vmatch; Jindx++){	//go through relevent J sequences

					length = const_d_numUniqueCharV[Vindx] + const_d_numUniqueCharJ[Jindx] + c_n;
					n_cnt = c_n;

					//check to see if sequence to create is the same length as the InVivo sequence
					if(seqLen != length) continue;

					////////////////////////////////////////////////////
					//begin generating sequences with no full chewbacks
					////////////////////////////////////////////////////

						sh_index = 3;											//reset our shared memory index
						seqMatch = 1;

						k = const_d_V_base[Vindx];									//starting address of V sequence
						for(int m = 0; m < const_d_numUniqueCharV[Vindx]; m++){		//go through each character in current V sequence
							tmpChar = const_d_V[k];									//load a V character into a temp variable
							if(tmpChar != iterSeq_sm[sh_index]){Vmatch = 0; break;}	//End V comparisons
							sh_index++;												//increment shared memory index
							k++;													//increment for next character
						}

						if(Vmatch == 0) continue;
				
						if(c_n != 0){
							//add n combination
							for(int m = 0; m < n_cnt; m++){
								tmpChar = nSeq[m]; 
								if(tmpChar != iterSeq_sm[sh_index]){seqMatch = 0; break;}	//jump to next iteration if a character does not match
								sh_index++;
							}
						}

						if(seqMatch == 0) continue;

						//glue current J sequence
						k = const_d_J_base[Jindx];								//starting address of V sequence
						for(int m = 0; m < const_d_numUniqueCharJ[Jindx]; m++){	//go through each character in current V sequence
							tmpChar = const_d_J[k];
							if(tmpChar != iterSeq_sm[sh_index]){seqMatch = 0; break;}	//jump to next iteration if a character does not match
							sh_index++;											//increment shared memory index
							k++;												//increment for next character
						}

						if(seqMatch == 0) continue;

						sum += c_DB_Full_Chew_Occur; //if we've made it this far, the sequences match.

						//Store Chewed back V Sequences
						//char * strcpy (d_chewArrV[d_counterIndex], const_d_V); 
						for (int i=0; i<const_d_numUniqueCharV[Vindx]; i++){
							d_chewArrV[d_counterIndex * 30 + i] = const_d_V[i];
						}

						//Store Chewed back J Sequences
						//char * strcpy (d_chewArrJ[d_counterIndex], const_d_J);
						for (int i=0; i<const_d_numUniqueCharJ[Jindx]; i++){
							d_chewArrJ[d_counterIndex * 30 + i] = const_d_J[i];
						}

						atomicAdd(&d_counterIndex, 1);  //Add to row iterator
						__syncthreads();

				} //end iterating through j sequences
		} //end iterating through v sequences



		//---------------------------------------------------------------------------------
		//Compare our InVivo Sequences to nDn combinations with V and J fully chewed
		//---------------------------------------------------------------------------------
		for(int Dindx = 0; Dindx < const_numDB1; Dindx++){			//go through all D sequences

			length = const_d_numUniqueCharDB1[Dindx] + c_n;
			n_cnt = c_n;
					
			//check to see if sequence to create is the same length as the InVivo sequence
			if(seqLen != length) continue;


			////////////////////////////////////////////////////
			//begin generating sequences with no full chewbacks
			////////////////////////////////////////////////////
			for(int j = 0; j < n_p1; j++){								//go through each n addition (n + 1)

				sh_index = 3;											//reset our shared memory index
				
				if(c_n != 0){
					//add n combination
					for(int m = 0; m < n_cnt; m++){
						tmpChar = nSeq[m]; 
						if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDn;}	//jump to next iteration if a character does not match
						sh_index++;
					}
				}

				//glue current D sequence
				k = const_d_DB1_base[Dindx];							//starting address of D sequence
				for(int m = 0; m < const_d_numUniqueCharDB1[Dindx]; m++){	//go through each character in current V sequence
					tmpChar = const_d_DB1[k];							//store V character in shared memory
					if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDn;}	//jump to next iteration if a character does not match
					sh_index++;											//increment shared memory index
					k++;												//increment for next character
				}

				if(c_n != 0){
					//add n combination
					for(int m = n_cnt; m < c_n; m++){
						tmpChar = nSeq[m]; 
						if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDn;}	//jump to next iteration if a character does not match
						sh_index++;
					}
				}

				n_cnt--; 

			
				sum += const_d_numOccurrenceDB1[Dindx]; //if we've made it this far, the sequences match.

				nDn: continue; //if there is no match go to next n variance
	    	} //end iterating through n sequences
	} //end iterating through d sequences




	//-----------------------------------------------------------------------------------
	//Compare our InVivo Sequences to nJ combinations with V and D chewed
	//-----------------------------------------------------------------------------------
	for(int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++){	//go through relevent J sequences

		length = const_d_numUniqueCharJ[Jindx] + c_n;
		n_cnt = c_n;
					

		//check to see if sequence to create is the same length as the InVivo sequence
		if(seqLen != length) continue;


		////////////////////////////////////////////////////
		//begin generating sequences with no full chewbacks
		////////////////////////////////////////////////////

		sh_index = 3;											//reset our shared memory index
				
		if(c_n != 0){
			//add n combination
			for(int m = 0; m < n_cnt; m++){
				tmpChar = nSeq[m]; 
				if(tmpChar != iterSeq_sm[sh_index]) goto nJ;	//jump to next iteration if a character does not match
					sh_index++;
				}
		}


		//glue current J sequence
		k = const_d_J_base[Jindx];								//starting address of V sequence
		for(int m = 0; m < const_d_numUniqueCharJ[Jindx]; m++){	//go through each character in current V sequence
			tmpChar = const_d_J[k];
			if(tmpChar != iterSeq_sm[sh_index]) goto nJ;	//jump to next iteration if a character does not match
			sh_index++;											//increment shared memory index
			k++;												//increment for next character
		}

		//sum++; //if we've made it this far, the sequences match. 
		sum += c_DB_Full_Chew_Occur;

		nJ: continue; //if there is no match go to next n variance

	} //end iterating through j sequences




	//----------------------------------------------------------------------------------
	//First compare our InVivo Sequences to nDnJ combinations with no full chewbacks
	//----------------------------------------------------------------------------------
		for(int Dindx = 0; Dindx < const_numDB1; Dindx++){			//go through all D sequences
			for(int Jindx = c_J_Begin; Jindx < c_J_End; Jindx++){	//go through relevent J sequences

				length = const_d_numUniqueCharJ[Jindx] + const_d_numUniqueCharDB1[Dindx] + c_n;
				n_cnt = c_n;
					

				//check to see if sequence to create is the same length as the InVivo sequence
				if(seqLen != length) continue;


				////////////////////////////////////////////////////
				//begin generating sequences with no full chewbacks
				////////////////////////////////////////////////////
				for(int j = 0; j < n_p1; j++){								//go through each n addition (n + 1)

					sh_index = 3;											//reset our shared memory index
				
					if(c_n != 0){
						//add n combination
						for(int m = 0; m < n_cnt; m++){
							tmpChar = nSeq[m]; 
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDnJ;}	//jump to next iteration if a character does not match
							sh_index++;
						}
					}

					//glue current D sequence
					k = const_d_DB1_base[Dindx];							//starting address of D sequence
					for(int m = 0; m < const_d_numUniqueCharDB1[Dindx]; m++){	//go through each character in current V sequence
						tmpChar = const_d_DB1[k];							//store V character in shared memory
						if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDnJ;}	//jump to next iteration if a character does not match
						sh_index++;											//increment shared memory index
						k++;												//increment for next character
					}

					if(c_n != 0){
						//add n combination
						for(int m = n_cnt; m < c_n; m++){
							tmpChar = nSeq[m]; 
							if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDnJ;}	//jump to next iteration if a character does not match
							sh_index++;
						}
					}

					//glue current J sequence
					k = const_d_J_base[Jindx];								//starting address of V sequence
					for(int m = 0; m < const_d_numUniqueCharJ[Jindx]; m++){	//go through each character in current V sequence
						tmpChar = const_d_J[k];
						if(tmpChar != iterSeq_sm[sh_index]){n_cnt--; goto nDnJ;}	//jump to next iteration if a character does not match
						sh_index++;											//increment shared memory index
						k++;												//increment for next character
					}

					n_cnt--; 

					//sum++; //if we've made it this far, the sequences match. 
					sum += const_d_numOccurrenceDB1[Dindx];

					nDnJ: continue; //if there is no match go to next n variance
	    		} //end iterating through n sequences
			} //end iterating through j sequences
		} //end iterating through d sequences


		//---------------------------------------------------------------------------------
		//Compare our InVivo Sequences to n with all full chewbacks
		//---------------------------------------------------------------------------------

		//check to see if sequence to create is the same length as the InVivo sequence
		if(seqLen != c_n) goto n_only_done;

		////////////////////////////////////////////////////
		//begin generating sequences with just n
		////////////////////////////////////////////////////
		sh_index = 3;											//reset our shared memory index

		//add n combination
		for(int m = 0; m < c_n; m++){
			tmpChar = nSeq[m]; 
			if(tmpChar != iterSeq_sm[sh_index]) goto n_only_done;	//jump to next iteration if a character does not match
			sh_index++;
		}

		//sum++; //if we've made it this far, the sequences match. 
		sum += c_DB_Full_Chew_Occur;

		n_only_done: 



		//-------------------------------------------------------------------------------------------------------
		//If only 1 thread-block, then we can write results to RAM using just InVivo sequence number
		//-------------------------------------------------------------------------------------------------------
		if(blockDim.x == 1){ //if there is only 1 thread per block, just use i as global memory index. No need for reduction
			d_Results[i] = sum;
		}


		//-------------------------------------------------------------------------------------------------------
		//					Perform Reduction of Results if more than 1 thread-block				  
		//-------------------------------------------------------------------------------------------------------
		//reduction for current InVivo sequence in shared memory
		if(blockDim.x > 1){
	
			result_sm[threadIdx.x] = sum;			//write a threads sum to the shared memory
			__syncthreads();						//make sure all sums have been written before proceeding

			int half = blockDim.x/2;

			while(1){								//how many reductions we need
				if(threadIdx.x < half){				//only certain threads perform reduction
					result_sm[threadIdx.x] += result_sm[threadIdx.x + half];
				}
				__syncthreads();
				if(half == 1) break;
				half = half/2;
			}

			__syncthreads();
			//write results to the global memory. Each thread-block writes 1 result for each InVivo Sequence i
			if(threadIdx.x == 0){					//we need only 1 thread in the thread block to write its result
				d_Results[i*gridDim.x + blockIdx.x] = result_sm[0];		//write our consolidated result into the global memory
			}

		} //end result reduction


	} //end iterating through InVivo Sequences
	return;
} //kernel done




#endif // #ifndef _TNT_KERNEL_H_
