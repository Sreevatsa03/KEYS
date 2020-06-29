//char pointers specified with _cp, integer pointers specified with _ip

/*================================================================================================*/
/*=											Definitions											 =*/
/*================================================================================================*/
#define N 10							//size of n in characters
#define N_BEGIN	0					//the smallest n value to test
#define NUM_V_FILES 20				//total number of V files
#define NUM_J_FILES 12				//total number of J files

#define NUM_CHAR_V_FILES 913		//total number of characters in all non-palandromic V files
#define NUM_CHAR_J_FILES 908		//total number of characters in all non-palandromic J files
//#define NUM_CHAR_V_FILES 3107		//total number of characters in all non-palandromic V files
//#define NUM_CHAR_J_FILES 3210		//total number of characters in all non-palandromic J files
#define NUM_V_SEQS 342				//total sequences in all non-palandromic V files
#define NUM_J_SEQS 271				//total sequences in all non-palandromic J files



/*================================================================================================*/
/*=										General Variables										 =*/
/*================================================================================================*/
char* gene_loc_cp;					//will contain the location of a specific gene file
int num_V_files = NUM_V_FILES;		//total number of V files
int num_J_files = NUM_J_FILES;		//total number of J files
int n   = 10;						//size of n
int threads;						//number of threads in a block
int blocks;							//number of thread-blocks per grid
unsigned int total_threads;			//total number of threads in a grid

/*================================================================================================*/
/*=										   DB1 Variables										 =*/
/*================================================================================================*/
//#define	memSizeDB1_1 1448	//number of characters for DB
#define	memSizeDB1_1 425	//number of characters for DB
#define	memSizeDB1_2 676	//169 DB sequences * sizeof(int)

int h_numCharFullDB1;				//number of characters in full chewback (should be 1), represented by X
char h_D1;							//First character in DB1 file, should be an X
int h_D1Occur;						//Number of occurrences for full chewback of DB1
int DB1size;						//number of characters in all DB1 chewbacks
int numDB1Unique;					//total number of unique DB1 sequences
int memSizeDB1;						//Size of memory to allocate for D chewbacks
int memSizeUniqueDB1;				//Size of memory to allocate for "numbers" of unique D sequences
char *h_DB1_cp;						//1D Array to hold DB1 chewbacks
int  *h_numOccurrenceDB1_ip;		//Number of times each D chewback occurs
int  *h_numUniqueCharDB1_ip;		//Number of characters in unique DB1 sequence occurence
int  *h_DB1_base_ip;				//gives starting index of D in h_DB1

/*================================================================================================*/
/*=										   InVivo Variables										 =*/
/*================================================================================================*/
int h_num_InVivo;					//total number of InVivo sequences
unsigned char *h_InVivo_cp64;		//contains InVivo sequences and other data with padding of 64
int InVivo_memSize64;				//size of memory required for sequences with padding length of 64

int *VJ_Pairs_ip = (int *)malloc(NUM_V_FILES * NUM_J_FILES * sizeof(int));	//contain number of sequences per VJ pair
int *VJ_Pair_Base_ip = (int *)malloc(NUM_V_FILES * NUM_J_FILES * sizeof(int));	//contain number of sequences per VJ pair
int VJ_Largest;						//largest VJ pair number

/*================================================================================================*/
/*=										 V Sequence Variables									 =*/
/*================================================================================================*/
char *h_V_cp = (char *)malloc(NUM_CHAR_V_FILES * sizeof(char));			//will contain sequences for all V files
int  *h_numUniqueCharV_ip = (int *) malloc(NUM_V_SEQS*sizeof(int));		//total number of characters in each V sequence for all V sequences in all V files
int  *h_V_base_ip = (int *)malloc(NUM_V_SEQS*sizeof(int));				//base address of each sequence of V
int  *numVUnique_ip = (int *)malloc(NUM_V_FILES*sizeof(int));			//total number of unique V Chewbacks for a given V file

int memSizeV = NUM_CHAR_V_FILES * sizeof(char);							//Total memory required to hold sequences in all V files
int memSizeOffsetV = NUM_V_SEQS*sizeof(int);							//total memory required to hold offset information for V files

/*================================================================================================*/
/*=										 J Sequence Variables									 =*/
/*================================================================================================*/
char *h_J_cp = (char *)malloc(NUM_CHAR_J_FILES * sizeof(char));			//will contain sequences for all J files
int  *h_numUniqueCharJ_ip = (int *) malloc(NUM_J_SEQS*sizeof(int));		//total number of characters in each J sequence for all J sequences in all J files
int  *h_J_base_ip = (int *)malloc(NUM_J_SEQS*sizeof(int));				//base address of each sequence of J
int  *numJUnique_ip = (int *)malloc(NUM_J_FILES*sizeof(int));			//total number of unique J Chewbacks for a given J file

int memSizeJ = NUM_CHAR_J_FILES * sizeof(char);							//Total memory required to hold sequences in all J files
int memSizeOffsetJ = NUM_J_SEQS*sizeof(int);							//total memory required to hold offset information for J files


