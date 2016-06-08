/*************************************************
misc.c
P. Clote
Miscellaneous routines for reading in RNA either as
string from command line or FASTA file, etc.

NOTES:
1) I chose a[] to be an array of short int, rather than
   char, since a[0] is the length of the original RNA
   sequence, stored in a[1],...,a[ a[0] ]
   Thus this requires a special irintRNA() function.

*************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>     //for INT_MAX
#include "misc.h"
#include <ctype.h>   	// character handling, eg toupper()
#include <stdbool.h>

void *xcalloc(size_t n, size_t s) {
  void *out = calloc(n,s);
  if (!out) {
    perror("calloc error");
    exit(1);
  }
  return out;
}

int auPair(char x, char y){
   return ( (x=='A'&&y=='U') || (x=='U'&&y=='A') || 
	    (x=='A'&&y=='T') || (x=='T'&&y=='A') );
}

int gcPair(char x, char y){
   return ( (x=='C'&&y=='G') || (x=='G'&&y=='C') );
}

int guPair(char x, char y){
   return ( (x=='G'&&y=='U') || (x=='U'&&y=='G') ||
	    (x=='G'&&y=='T') || (x=='T'&&y=='G') );
}

int watsonCrick(char x, char y){
   return( auPair(x,y) || gcPair(x,y) );
}

int basePair(char x, char y){
   return( watsonCrick(x,y) || guPair(x,y) );
}

void printSecStrListAsString(char *secStr){
  /* Print the secondary structure as a string */
  int i=0;
  while (secStr[i]!='\0') printf("%c",secStr[i++]);
  printf("\n");
}

void printBasePairList(int *bps, int len){
  /* Print the secondary structure as a string */
  int i;
  for (i=0;i<len;i++) printf("(%u,%u) ",bps[2*i],bps[2*i+1]);
  printf("\n");
}

int energyOfBasePair(char x, char y) {
  if (auPair(x,y))   return AUenergy; 
  else if (guPair(x,y)) return GUenergy;
  else if (gcPair(x,y)) return GCenergy;
  else return 0;
}

int *getBasePairList(char *secStr) {
  /* Returns list L of ordered pairs (i,j) where i<j and
   * positions i,j occupied by balancing parentheses
   * For linear time efficiency, use stack
   * Assume that secStr is string consisting of '(',')' and '.'
   * Values -2,-1 returned mean NOT well balanced
   * -2 means too many ( with respect to )
   * -1 means too many ) with respect to (
   * If 1,-1 not returned, then return (possibly empty) list */
  
  int len = strlen(secStr);
  int *S = (int *) xcalloc(len/2,sizeof(int));  //empty stack
  int *L = (int *) xcalloc(2*len*(len-1)/2+1, sizeof(int)); /* initially empty
							     * list of base 
							     * pairs */
  int j, k = 0, l = 1;
  char ch;

  for (j=0;j<len;j++) {
    ch = secStr[j];
    if (ch == '(')
      S[k++] = j;
    else if (ch == ')') {
      if (k==0) {
	L[0] = -1; 
	return L;
      }
      else {
        L[2*(l)] = S[--k];
	L[2*(l++)+1] = j;
      }
    }
  }

  /* Store the number of base pairs in the array */
  L[2*l] = -1;

  if (k != 0) {
    L[0] = -2;
  }
  
  free(S);

  return L;
}

int validNucleotideString(char *rna) {
  int i=-1;
  while (rna[++i]!='\0') {
    rna[i] = toupper(rna[i]);
    if (rna[i]=='T')
      rna[i]='U';
    else if (rna[i]!='A' && rna[i]!='G' && rna[i]!='C' && rna[i]!='U')
      return 0;
    }
  return 1;
}

int countmatches(char *rna, char ch) {
  int i=-1, count = 0;
  while (rna[++i]!='\0') {
    if (rna[i]==ch)
      count++;
  }
  return count;
}

int numParenPairs(char *secStr) {
  /* Assume that secStr is well-balanced parenthesis string with dots */
  return  countmatches(secStr,'(');
}

char *getInputRNA( char *fileName){
	/* Function to retrieve RNA sequence from a fasta file
	* or as FASTA file.
	* WARNING: rna is 0-indexed.  */
	FILE * infile;
	char ch;
	int i,n,len;
	char *rna, tmprna[NN]; 

	/* Read FASTA file of one RNA sequence  */ 
	infile = fopen(fileName,"r");
	ch=fgetc(infile);
	if (ch == '>')
		while (ch != '\n') ch = fgetc(infile);
	n=0; //current position of insert into RNA string
	while (ch != EOF && n < NN){
		if (ch == '\n') {
		ch = fgetc(infile);
		if (ch == EOF) { tmprna[n]='\0'; break;	}
	}
		tmprna[n] = toupper(ch);
		n++;
		ch = fgetc(infile); 
	}
	tmprna[n] = '\0'; //string terminator
	if (n>=NN){
	  fprintf(stderr,"Length of RNA exceeds array size %d\n",NN);
	  exit(1);
	}
	fclose(infile);
	rna = (char *) xcalloc(n+1,sizeof(char));
	n=0;
	while (tmprna[n]!='\0') {
		rna[n] = tmprna[n];
		n++;
	}
	rna[n]='\0';
	return rna;
} 

/* Useful functions for statistics */

int max(int *l, int n) {
  int m = 0, i;
  for (i=0;i<n;i++)
    if (m<l[i])
      m = l[i];
  return m;
}

int min(int *l, int n) {
  int m = INF, i;
  for (i=0;i<n;i++)
    if (m>l[i])
      m = l[i];
  return m;
}

double mean(int *l, int n) {
  int s = 0, i;
  for (i=0;i<n;i++)
    s += l[i];
  return ((double) s)/n;
}

double std(int *l, int n) {
  double s = 0.0;
  int i;
  double m = mean(l,n);
  for (i=0;i<n;i++)
    s += (l[i] - m)*(l[i] - m);
  return sqrt(s/n);
}

double dmean(double *l, int n) {
  double s = 0.0;
  int i;
  for (i=0;i<n;i++)
    s += l[i];
  return s/n;
}

double dstd(double *l, int n) {
  double s = 0.0;
  int i;
  double m = dmean(l,n);
  for (i=0;i<n;i++)
    s += (l[i] - m)*(l[i] - m);
  return sqrt(s/n);
}

/* Random integer in interval a..b */
int randint(int a, int b){
  int n, r;
  if (b<a)
    return 0;
  n = b-a+1;
  r = lrand48();
  if (n==0)
    return a;
  else
    return r % n + a;
}

/* Base pairing matrix */
bool *canBasePair(char *rna) {
  /* WARNING: rna is 0-indexed. 
   * Returns BP[i*n+j] = 0 or 1. */
  int i, j, n = strlen(rna);
  bool *BP = (bool *) xcalloc((n+1)*(n+1),sizeof(bool));
  for (i=0;i<n-1-THRESHOLD;i++) {
    for (j=i+THRESHOLD+1;j<n;j++)
      if (basePair(rna[i],rna[j]))
	BP[i*n+j] = 1;
  }
  return BP;
}

/* Number of base pairs in base pair list (the list must be ended with -1) */
int length(int *L) {
  int n=0;
  while (L[2*n]!=-1) {
    n++;
  }
  return n;
}

char* getExecPath(char* argv0){
	char* exec_path=(char*) malloc (sizeof(char)*PATH_MAX);
	exec_path[0]='.';
	exec_path[1]='\0';
	int path_length=1;
	FILE* file;
	#if defined(_WIN32)
	// GetModuleFileName
	//_pgmptr
	#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	char buff[PATH_MAX];
	int bufsize = PATH_MAX-1;
	if(file = fopen("/proc/self/exe", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/exe", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/curproc/file", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/curproc/file", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/self/path/a.out", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/path/a.out", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else{
		exec_path= argv0;
	}
	char* slash_pointer= strrchr(exec_path,'/');
	int slash_position = (int)(slash_pointer - exec_path);
	if(slash_position != path_length){			
		exec_path[slash_position]='\0';
	}
	else{
		exec_path[0]='.';
		exec_path[1]='\0';
	}

	#endif
	

	return exec_path;
}
