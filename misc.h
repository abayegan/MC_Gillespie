#ifndef MISC_H
#define MISC_H

#define DEBUG 0

/* Constant definitions */
#define THRESHOLD  3
#define INF INT_MAX/2
#define AUenergy  -2
#define GUenergy  -1
#define GCenergy  -3
#define NN 610 /* Maximum length of RNA */

#include <stdbool.h>

void *xcalloc(size_t , size_t );
int guPair(char, char);
int auPair(char, char);
int gcPair(char, char);
int watsonCrick(char, char);
int basePair(char, char);
void printSecStrListAsString(char *);
void printBasePairList(int *, int);
int *getBasePairList(char *);
int validNucleotideString(char *);
int countmatches(char *, char);
int numParenPairs(char *);
char *getInputRNA(char *);
double max(double *, int);
double min(double *, int);
double mean(double *, int);
double std(double *, int);
double dmean(double *, int);
double dstd(double *, int);
int randint(int,int);
/* Base pairing matrix */
bool *canBasePair(char *);
int length(int *);
void printUsage(char * );
void printHelp(char *);
char* getExecPath(char*);
#endif
