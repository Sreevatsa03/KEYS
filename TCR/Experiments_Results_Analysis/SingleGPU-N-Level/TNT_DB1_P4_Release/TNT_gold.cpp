#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include <ios>
#include <iostream>
#include <cstdint>
#include <memory.h>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
// export C interface
extern "C"
void msgs( int);

extern "C" //get file location of DB1
char *getDB1FileLoc();

extern "C" //retrieve DB1 sequences from file
void getDB1( int*, int*, int*, int*, char**, int**, int**, int**, char*, int*, int*);

extern "C" //print DB1 data to a test file
void printD( int*, int*, int*, int*, char**, int**, int**, int**, char*, int*, int*);

extern "C" //retrieve V sequences from files
void getV(char*, int*, int*, int*);

extern "C" //retrieve J sequences from files
void getJ(char*, int*, int*, int*);

extern "C" //get file location of a V file
char *getVFileLoc( int);

extern "C" //get file location of a J file
char *getJFileLoc( int);

extern "C" //print V sequences to a file, to check correctness
void printV(char*, int*);

extern "C" //print J sequences to a file, to check correctness
void printJ(char*, int*);

extern "C" //Get InVivo sequence fro file
void getInVivo( int*, int*, unsigned char**);

extern "C" //Get InVivo sequence fro file
void printInVivo( int*, int*, unsigned char**);

extern "C" //Get number of pairs in each VJ combinations
void getNum_VJ_Pairs( int, int, int*, int*);

extern "C" //print the InVivo Results
void print_InVIvo_Results(int, int, int, int, unsigned int*, int);

////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////
//Messages
/////////////////////////////////////////////////
void
msgs(int i){

	switch(i){
			case 0: printf("\nERROR: This size of 'n' is not supported.\n'n' must be between 1 and 9\n\n"); break;
			case 1: printf("\nERROR: DB1 File Could Not be Opened \n"); break;
			case 2: printf("\nDB1 File Opened Successfully\nRetrieving DB1 Information\n\n"); break;
			case 3: printf("\nERROR: No case of fully chewed back DB1\n\n"); break;
			case 4: printf("\nERROR: File number does not match for V_Genes\n\n"); break;
			case 5: printf("\nERROR: File number does not match for J_Genes\n\n"); break;
			case 6: printf("\nERROR: V File could not be opened\n\n"); break;
			case 7: printf("\nERROR: V Palandromic output file could not be opened\n\n"); break;
			case 8: printf("\nERROR: V non-Palandromic output file could not be opened\n\n"); break;
			case 9: printf("\nERROR: J File could not be opened\n\n"); break;
			case 10: printf("\nERROR: J Palandromic output file could not be opened\n\n"); break;
			case 11: printf("\nERROR: J non-Palandromic output file could not be opened\n\n"); break;
	}

}

/////////////////////////////////////////////////
//returns the file location for a DB1 sequences
/////////////////////////////////////////////////
char*
getDB1FileLoc(){
	char* gene_loc;

	gene_loc = "Data/D_Genes/DB1_Chewback_Info.txt";

	return gene_loc;
}


/////////////////////////////////////////////////
//returns the file location for a V sequences
/////////////////////////////////////////////////
char*
getVFileLoc(int i){

	if(i < 0 || i > 19){msgs(4); exit(1);}

	char* gene_loc;

	switch(i){
		   case 0:gene_loc = "Data/V_Genes/1.txt"; break;
		   case 1:gene_loc = "Data/V_Genes/2.txt"; break;
		   case 2:gene_loc = "Data/V_Genes/3.txt"; break;
		   case 3:gene_loc = "Data/V_Genes/4.txt"; break;
		   case 4:gene_loc = "Data/V_Genes/5.txt"; break;
		   case 5:gene_loc = "Data/V_Genes/6.txt"; break;
		   case 6:gene_loc = "Data/V_Genes/7.txt"; break;
		   case 7:gene_loc = "Data/V_Genes/8.txt"; break;
		   case 8:gene_loc = "Data/V_Genes/9.txt"; break;
		   case 9:gene_loc = "Data/V_Genes/10.txt"; break;
		   case 10:gene_loc = "Data/V_Genes/11.txt"; break;
		   case 11:gene_loc = "Data/V_Genes/12.txt"; break;
		   case 12:gene_loc = "Data/V_Genes/13.txt"; break;
		   case 13:gene_loc = "Data/V_Genes/14.txt"; break;
		   case 14:gene_loc = "Data/V_Genes/15.txt"; break;
		   case 15:gene_loc = "Data/V_Genes/16.txt"; break;
		   case 16:gene_loc = "Data/V_Genes/17.txt"; break;
		   case 17:gene_loc = "Data/V_Genes/18.txt"; break;
		   case 18:gene_loc = "Data/V_Genes/19.txt"; break;
		   case 19:gene_loc = "Data/V_Genes/20.txt"; break;
	}

	return gene_loc;
}


/////////////////////////////////////////////////
//returns the file location for a J sequences
/////////////////////////////////////////////////
char*
getJFileLoc(int i){

	if(i < 0 || i > 11){msgs(5); exit(1);}

	char* gene_loc;

	switch(i){
		   case 0:gene_loc = "Data/J_Genes/1.txt"; break;
		   case 1:gene_loc = "Data/J_Genes/2.txt"; break;
		   case 2:gene_loc = "Data/J_Genes/3.txt"; break;
		   case 3:gene_loc = "Data/J_Genes/4.txt"; break;
		   case 4:gene_loc = "Data/J_Genes/5.txt"; break;
		   case 5:gene_loc = "Data/J_Genes/6.txt"; break;
		   case 6:gene_loc = "Data/J_Genes/7.txt"; break;
		   case 7:gene_loc = "Data/J_Genes/8.txt"; break;
		   case 8:gene_loc = "Data/J_Genes/9.txt"; break;
		   case 9:gene_loc = "Data/J_Genes/10.txt"; break;
		   case 10:gene_loc = "Data/J_Genes/11.txt"; break;
		   case 11:gene_loc = "Data/J_Genes/12.txt"; break;
	}

	return gene_loc;
}


/////////////////////////////////////////////////
//gets all DB1 sequences
/////////////////////////////////////////////////
void
getDB1( int* numDB1Unique, int* DB1size, int* memSizeDB1, int* memSizeUniqueDB1,
        char** h_DB1_cp, int** h_numOccurrenceDB1_ip, int** h_numUniqueCharDB1_ip, int** h_DB1_base_ip,
		char* gene_loc_cp, int* h_numCharFullDB1, int* h_D1Occur){

	char h_D1;					//First character in DB1 file, should be an X
	int i, j, k, sum = 0;		//incrementers

	//create a pointer to the file containing DB1
	fstream IN_FILE;
	IN_FILE.open (gene_loc_cp, ios::in);

	//Print error message if file cannot be opened
	if(IN_FILE == NULL){msgs(1); exit(1);}

	//Print success message if file was opened
	if(IN_FILE != NULL)msgs(2);

	IN_FILE >> *numDB1Unique;		//Retrieve total number of unique DB1 chewbacks
	IN_FILE >> *DB1size;			//Retrieve total size of DB1 chewbacks in terms of characters

	*DB1size -= 1;					//Remove 1 character for full chewback
	*numDB1Unique -= 1;				//Remove 1 to account for full chewback

	// Determine size of memory to allocate for D chewbacks
	*memSizeDB1 = *DB1size * sizeof(char);

    // Determine size of memory to allocate for "numbers" of unique D sequences
    *memSizeUniqueDB1 = *numDB1Unique * sizeof(int);

	//Allocate memory on host for 1D Array Containing DB1 chewbacks and number of chewbacks
	// and number of characters in a unique occurence of DB1, and gives starting index of D in h_DB1
    *h_DB1_cp              = (char* )malloc(*memSizeDB1);
    *h_numOccurrenceDB1_ip = (int * )malloc(*memSizeUniqueDB1);
    *h_numUniqueCharDB1_ip = (int * )malloc(*memSizeUniqueDB1);
    *h_DB1_base_ip         = (int * )malloc(*memSizeUniqueDB1);

	 //assume first chewback is 'X' or fully chewed back
	IN_FILE >> *h_numCharFullDB1;		//number of characters in full chewback (should be 1)
	IN_FILE >> h_D1;					//get first character to make sure it is X
	if(h_D1 != 'X'){msgs(3); exit(1); } //if no full chewback report error
	IN_FILE >> *h_D1Occur;				//get the number of full chewback occurrences

	 //get DB1 Chewbacks and number of each chewback from file
	 //here we are getting all of the rest of information from the input file
	k = 0;
	for(i=0; i < *numDB1Unique; i++){						//increment through each of the unique sequences in the file
		IN_FILE >> (*h_numUniqueCharDB1_ip)[i];				//get number of characters in this sequence
		(*h_DB1_base_ip)[i] = sum;							//store running sum of each sequence length
		sum += (*h_numUniqueCharDB1_ip)[i];					//increment sum for determining base address of next sequence
		for(j = 0; j < (*h_numUniqueCharDB1_ip)[i]; j++){	//get the unique sequence, store in h_DB1
			IN_FILE >> (*h_DB1_cp)[k];						//store character from unique sequence in h_DB1
			k++;											//incrememnt k for next character in h_DB1
		}
		IN_FILE >> (*h_numOccurrenceDB1_ip)[i];				//get the number of unique number of occurrences of this chewback
	 }
 	IN_FILE.close();										//close file for DB1 on hard disk
	return;
}

/////////////////////////////////////////////////
//prints DB1 sequences to a file
/////////////////////////////////////////////////
void
printD( int* numDB1Unique, int* DB1size, int* memSizeDB1, int* memSizeUniqueDB1,
        char** h_DB1_cp, int** h_numOccurrenceDB1_ip, int** h_numUniqueCharDB1_ip, int** h_DB1_base_ip,
		char* gene_loc_cp, int* h_numCharFullDB1, int* h_D1Occur){

	int i, j, k;				//incrementers
	int total_SeqD;				//total number of D sequences in either Palandromic or non-Palandromic

	total_SeqD = *numDB1Unique;	//set total number of D sequences

	fstream OUT_FILE;			//create a pointer to the file to output D

	OUT_FILE.open("Result_out/Ddata_nonPalandromic.txt", ios::out);	//open file for output
	if(OUT_FILE == NULL){printf("Could not open D output file\n\n"); exit(1);}						//print error message if V file could not be opened


	k = 0;
	for(i = 0; i < total_SeqD; i++){								//iterate through each of the V sequences

		for(j = (*h_DB1_base_ip)[i]; j < (*h_numUniqueCharDB1_ip)[i] + (*h_DB1_base_ip)[i]; j++){				//iterate through number of characters in current D sequence
			OUT_FILE << (*h_DB1_cp)[k]; k++;								//write new V character to output file
		}

		OUT_FILE << "\n";											//create new line
	}

 	OUT_FILE.close();

	return;



}


/////////////////////////////////////////////////
//gets all V sequences from all V files
/////////////////////////////////////////////////
void
getV(char* h_V_cp, int* h_numUniqueCharV_ip, int* h_V_base_ip, int* numVUnique_ip){

	int i, j, k, l, m, sum;	//incrementers
	char* gene_loc;			//will contain the file location for each V file

	//create a pointer to the file containing V
	fstream IN_FILE;
	k = 0, l = 0, m = 0, sum = 0;
    //we will iterate through each of the 20 V files
	for(i = 0; i < 20; i++){
		gene_loc = getVFileLoc(i);							//determine current V file location
		IN_FILE.open (gene_loc, ios::in);					//open the file on disk
		if(IN_FILE == NULL){msgs(6); exit(1);}				//print error message if V file could not be opened
		IN_FILE >> numVUnique_ip[i];						//Read the number of unique V sequences in this file
		for(j = 0; j < numVUnique_ip[i]; j++){				//iterate through each of the V sequences in this file
			IN_FILE >> h_numUniqueCharV_ip[m];				//read in the total number of characters contained in the current sequence
			h_V_base_ip[m] = sum;							//store running sum of each sequence length
			sum += h_numUniqueCharV_ip[m];					//increment sum for next sequence base address
			for(k = 0; k < h_numUniqueCharV_ip[m]; k++){	//iterate through each of the characters in the current sequence
				IN_FILE >> h_V_cp[l];						//store sequence into array
				l++;
			}m++;
		}
 		IN_FILE.close();
	}
	return;
}



//////////////////////////////////////////////////////
//print V sequences to a file, to check correctness
//////////////////////////////////////////////////////
void
printV(char* h_V_cp, int* h_numUniqueCharV_ip){

	int i, j, k;			//incrementers
	int total_SeqV;			//total number of V sequences in either Palandromic or non-Palandromic

	fstream OUT_FILE;		//create a pointer to the file to output V


	OUT_FILE.open("Result_out/Vdata.txt", ios::out);	//open file for output
	total_SeqV = 342;											//set total number of V sequences
	if(OUT_FILE == NULL){msgs(8); exit(1);}						//print error message if V file could not be opened


	k = 0;
	for(i = 0; i < total_SeqV; i++){								//iterate through each of the V sequences
		for(j = 0; j < h_numUniqueCharV_ip[i]; j++){				//iterate through number of characters in current V sequence
			OUT_FILE << h_V_cp[k]; k++;								//write new V character to output file
		}
		OUT_FILE << "\n";											//create new line
	}

 	OUT_FILE.close();

	return;
}


/////////////////////////////////////////////////
//gets all J sequences from all J files
/////////////////////////////////////////////////
void
getJ(char* h_J_cp, int* h_numUniqueCharJ_ip, int* h_J_base_ip, int* numJUnique_ip){

	int i, j, k, l, m, sum;	//incrementers
	char* gene_loc;			//will contain the file location for each V file

	//create a pointer to the file containing V
	fstream IN_FILE;
	k = 0, l = 0, m = 0, sum = 0;
    //we will iterate through each of the 20 V files
	for(i = 0; i < 12; i++){
		gene_loc = getJFileLoc(i);							//determine current J file location
		IN_FILE.open (gene_loc, ios::in);					//open the file on disk
		if(IN_FILE == NULL){msgs(9); exit(1);}				//print error message if J file could not be opened
		IN_FILE >> numJUnique_ip[i];						//Read the number of unique J sequences in this file
		for(j = 0; j < numJUnique_ip[i]; j++){				//iterate through each of the J sequences in this file
			IN_FILE >> h_numUniqueCharJ_ip[m];				//read in the total number of characters contained in the current sequence
			h_J_base_ip[m] = sum;							//store running sum of each sequence length
			sum += h_numUniqueCharJ_ip[m];					//increment sum for next sequence base address
			for(k = 0; k < h_numUniqueCharJ_ip[m]; k++){	//iterate through each of the characters in the current sequence
				IN_FILE >> h_J_cp[l];						//store sequence into array
				l++;
			}m++;
		}
 		IN_FILE.close();
	}
	return;
}


//////////////////////////////////////////////////////
//print J sequences to a file, to check correctness
//////////////////////////////////////////////////////
void
printJ(char* h_J_cp, int* h_numUniqueCharJ_ip){

	int i, j, k;			//incrementers
	int total_SeqJ;			//total number of J sequences in either Palandromic or non-Palandromic

	fstream OUT_FILE;		//create a pointer to the file to output J


	OUT_FILE.open("Result_out/Jdata.txt", ios::out);	//open file for output
	total_SeqJ = 271;											//set total number of J sequences
	if(OUT_FILE == NULL){msgs(11); exit(1);}					//print error message if J file could not be opened


	k = 0;
	for(i = 0; i < total_SeqJ; i++){								//iterate through each of the J sequences
		for(j = 0; j < h_numUniqueCharJ_ip[i]; j++){				//iterate through number of characters in current J sequence
			OUT_FILE << h_J_cp[k]; k++;								//write new J character to output file
		}
		OUT_FILE << "\n";											//create new line
	}

 	OUT_FILE.close();

	return;
}



/////////////////////////////////////////////////
//gets all InVivo Sequences
/////////////////////////////////////////////////
void
getInVivo( int* InVivo_memSize64, int* h_num_InVivo, unsigned char** h_InVivo_cp64){

	int i, ii, j, k;
	int numChar;					//number of characters in the file
	int whichSeq;
	int temp;
	unsigned char temp2;

	//create a pointer to the file containing InVivo Sequences
	fstream IN_FILE;
	IN_FILE.open ("Data/InVivo/All_In_Vivo_Nucleotypes.txt", ios::in);

	//Print error message if file cannot be opened
	if(IN_FILE == NULL){msgs(12); exit(1);}


	IN_FILE >> *h_num_InVivo;		//Retrieve total number of InVivo Sequences
	IN_FILE >> numChar;				//Retrieve total number of sequence characters


	printf("num InVivo = %d\n", *h_num_InVivo);

	// Determine size of memory to allocate for InVivo Sequences with length 64, 128, and 256 including padding
	*InVivo_memSize64  = *h_num_InVivo *  64 * sizeof(char);

	//Allocate host memory for InVivo sequences and other data
    *h_InVivo_cp64  = (unsigned char* )malloc(*InVivo_memSize64);

	printf("\nTotal memory needed for InVivo Sequences \nwith three versions of padding: %1.2f MB\n\n",
		   (float)(*InVivo_memSize64)/ 1048576);


	//get each of the InVivo sequences and add padding
	for(i=0; i<*h_num_InVivo; i++){

		ii = 64*i;

		IN_FILE >> whichSeq;										//which sequence are we looking at
		IN_FILE >> temp; (*h_InVivo_cp64)[ii]  = temp; ii++;		//Read in V value

		IN_FILE >> temp; (*h_InVivo_cp64)[ii]  = temp; ii++;		//Read in J value

		IN_FILE >> temp; (*h_InVivo_cp64)[ii]  = temp; ii++;		//Read in Length

		for(j = 0; j < temp; j++){ //get all the characters for our current sequence
			IN_FILE >> temp2; (*h_InVivo_cp64)[ii]  = temp2; ii++;
		}

		//add padding up to a length of 64
		for(k = temp+3; k<64; k++){ //only want to add padding beyond sequence information
			(*h_InVivo_cp64)[ii]  = 0; ii++;
		}

	}

 	IN_FILE.close();										//close file for DB1 on hard disk

	return;
}

/////////////////////////////////////////////////
//prints all InVivo Sequences
/////////////////////////////////////////////////
void
printInVivo( int* InVivo_memSize64, int* h_num_InVivo, unsigned char** h_InVivo_cp64){


	int i, j, k, ii, jj, kk;				//incrementers
	int total_Seq;				//total number of sequences

	total_Seq = *h_num_InVivo;	//set total number of sequences

	fstream OUT_FILE64;			//create pointers to files

	OUT_FILE64.open("Result_out/InVivo64.txt",   ios::out);	//open files for output

	if(OUT_FILE64  == NULL){printf("Could not open InVivo64 output file\n\n");  exit(1);}						//print error message if V file could not be opened

	k = 0; ii = 0; kk = 0; jj = 0;
	for(i = 0; i < total_Seq; i++){						//iterate through each of the V sequences

		OUT_FILE64 << (int)(*h_InVivo_cp64)[ii] << " "; ii++;
		OUT_FILE64 << (int)(*h_InVivo_cp64)[ii] << " "; ii++;
		OUT_FILE64 << (int)(*h_InVivo_cp64)[ii] << " "; ii++;


		for(j = 3; j < 64; j++){				//iterate through number of characters in current D sequence
			OUT_FILE64 << (*h_InVivo_cp64)[ii]; ii++;			//write new V character to output file
		}


		OUT_FILE64  << "\n";
	}


 	OUT_FILE64.close();

	return;

}


/////////////////////////////////////////////////
//gets the number of sequences per VJ pair
/////////////////////////////////////////////////
void
getNum_VJ_Pairs( int NUM_V_FILES, int NUM_J_FILES, int* VJ_Pairs_ip, int* VJ_Largest){

	*VJ_Largest = 0;

	fstream IN_FILE;

	IN_FILE.open("Data/InVivo/VJ_Pairs.txt", ios::in);

	int j = NUM_V_FILES * NUM_J_FILES;

	for(int i=0; i< j; i++){
		IN_FILE >> VJ_Pairs_ip[i];
		if(VJ_Pairs_ip[i] > *VJ_Largest)
			*VJ_Largest = VJ_Pairs_ip[i];
	}

	IN_FILE.close();

	return;
}



////////////////////////////////////////////////////////////
//determine which kernel to run
////////////////////////////////////////////////////////////
extern "C" //print the InVivo Results
void print_InVIvo_Results(int n, int V, int J, int numResults, unsigned int* h_Result, int InVivoBase){

	int i, j;


	//local variables
	char buff_v[8];
	char buff_j[8];
	char buff_n[8];

	//convert input values to characters
	
	int n_v = sprintf(buff_v, "%d", ( V + 1 ));
	int n_j = sprintf(buff_j, "%d", ( J + 1 ));
	int n_n = sprintf(buff_n, "%d", n);

	
	string ResFile = "Result_out/InVivo/n_";
	ResFile.append(buff_n); 		//append n number
	ResFile.append("_v_");
	ResFile.append(buff_v); 		//append v number
	ResFile.append("_j_");
	ResFile.append(buff_j); 		//append j number
	ResFile.append(".txt");

	//convert string to character array for Key file
	char *ResFile_location=new char[ResFile.size()+1];
	ResFile_location[ResFile.size()]=0;
	memcpy(ResFile_location, ResFile.c_str(),ResFile.size());

	fstream OUT_FILE;								//create a pointer to the file to output

	OUT_FILE.open(ResFile_location, ios::out);		//open file for output

	for(i = 0, j = InVivoBase; i < numResults; i++, j++){
		OUT_FILE << j+1 << "\t" << h_Result[i] << "\n";			//write new V character to output file
	}

	OUT_FILE.close();		//close

	return;
}

