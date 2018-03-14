#include <stdio.h>
#include <stdlib.h>     //for RAND_MAX
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <limits.h>     //for INT_MAX
#include "misc.h"
#include <stdbool.h>


/* From libRNA.a */
double temperature;
int dangles;
float energy_of_struct(char*, char*);
float fold(const char *, char *);

/* Prototypes */
char *initialEmptySecStr(char *);
int canAddBasePair(int, int, int, int *);
void listOfNextStatesByAddingBasePair(char *, char *, int *,int **);
void listOfNextStatesByAddingCompatibleBasePair(char *, char *, int *,int **);
struct stateConfig monteCarlo(char *);
void nextConfig(char *, struct stateConfig *, int *);
int run(char *, int);
inline double sampleExpDist(double);
void printUsage(char * );
void printHelp(char *);


struct stateConfig {
	int num;
	double time;
	char *secStr;
	float MFE;
	int *bps;
	int folded;
};

#define MFPT 0

/* Sample from the exponential distribution with given mean */
inline double sampleExpDist(double mean){
  return -mean*log(drand48());
  }
/* Defined constants */
#define MIN(x,y) ((x)<(y)?x:y)
#define ABS(x) ((x)>=0?(x):-(x))
/* Global variables */
int trjLength=-1;  /* number of Monte Carlo steps, -1 indicates stop at the MFE regardless of the number of steps */
bool *BP;       // BP[i*n+j]=0 or 1 if (i,j) is potential base pair 
FILE *timeFile; // file var for output file of folding times
float minE,RT= 0.00198717 * 310.15;     // minimum free energy returned by Vienna fold() 
char turner99[100] = "/param_files/rna_turner1999.par";
char turner04[100] = "/param_files/rna_turner2004.par";
char andronescu07[100] = "/param_files/rna_andronescu2007.par";
char method[]="gil",engModel[] = "04"; // time-driven MC, event drive MC or Gillespie algorithms
int trjLength, seed, numRun=1, synchronized=0, verbose=0, hasting=0;
double MFE; /*minimum free energy of the input sequence*/
char * mfeSS; /*minimum free energy structure of the input sequence*/
char * targetStr = NULL; /*target structure to stop the trajectories*/
int GO = 0;
void printUsage(char *name){
	printf("\nUSAGE: %s [-s sequence]|[-f FastaFile] [-m method] [-y random seed] [-r number of runs] [-l trajectory length] [-d 0|1|2] [-t temperature(c)] [-e 99|04] [-o output file prefix ] [-v]\ntry %s -h for help\n\n",name, name);
	exit(1);
}

void printHelp(char *name){
	printf("\nUSAGE: %s [-s sequence]|[-f FastaFile] [-m method] [-y random seed] [-r number of runs] [-l trajectory length] [-d 0|1|2] [-t temperature(c)] [-e 99|04] [-o output file prefix ] [-v]\ntry %s -h for help\n\n",name, name);
	printf("This program simulates kinetics of folding with event-driven Monter Carlo or Gillespie algorithms. The available options are as follows:\n"\
	"-s : input RNA sequence\n"\
	"-f : input fasta file containing the RNA. This flag and '-s' are mutually exclusive\n"\
	"-m : Define which method should be used to produce trajectories. Use \"emc\" or \"gil\" for event-driven monte carlo or gillespie respecively\n"\
	"-y : fix the random seed to the given value. This allows to compare the results from the methods\n"\
	"-r : define the number of trajectories for the given RNA. The default is 1.\n"\
	"-l : define the number of states after which trajectory stops. The default is stopping at the MFE structure\n"\
	"-d : set the dangles to be 0,1 or 2 to compute the energy of each structures\n"\
	"-t : set the temprature (C) to compute the energy of structures\n"\
	"-e : set the energy model to be Turner 99, Turner 2004 or Andronescu 2007  by entering 99|04|07 respectively\n"\
	"-o : set the output folder prefix name to save the mean first passage times. The default prefix is 'rna'. The prefix will be followed by '_firstPassageTimes'\n"\
	"-v : print the trajectory on the output. Save the printed output if you need the trajectories. The output fields are (1)the state,(2)energy of state, (3)time stayed in the state, (4)number of neighbors of the state\n");
	exit(1);
}

/* 
 * Create an initial secondary structure of the same length as rna
 * with no base pairs. Both rna and secStr are 0-indexed. 
 */
char *initialEmptySecStr(char *rna) {
  int i=-1, n = strlen(rna);
  char *secStr = (char *) xcalloc(n+1,sizeof(char));
  while (rna[++i]!='\0')
    secStr[i] = '.';
  return secStr;
}

int canAddBasePair(int i, int j, int n, int *bps) {
  /* n is length of RNA sequence with secondary structure secStr.
   * bps is list of base pairs of secStr obtained by function
   * getBasePairList() or obtained incrementally from an empty list
   * Current function tests if 0<= i < j < n and whether for all (x,y)
   * in bps it is NOT the case that either i<x<j<y or x<i<y<j.  
   * WARNING:
   * This function does NOT check that nucleotides in positions i,j
   * can actually base pair. This is done by basePair(rna[i],rna[j])
   */
  int k=0, x, y;
  assert ((0<=i)&&(i<j)&&(j<n));
  while (bps[2*k]!=-1) {
    x = bps[2*k]; y = bps[2*k+1];
    if ( ((i<x) && (x<j) && (j<y)) ||
	 ((x<i) && (i<y) && (y<j)) ) 
      return 0;
    k++;
  }
  return 1;
}

  /* bps is list of base pairs from secStr, a secondary structure for rna
   * For base pair (i,j) that can be added, will assign:
   *          *basePairsWhichCanBeAdded[2*k]=i
   *          *basePairsWhichCanBeAdded[2*k+1]=j
   * WARNING: Base pairs are 0-indexed.
   *          *basePairsWhichCanBeAdded[2*k]= -1
   * indicates that the last base pair has been reached. */
void listOfNextStatesByAddingBasePair(char *rna, char *secStr, int *bps, 
				      int **basePairsWhichCanBeAdded) {
  int i,j,k = 0;
  int n = strlen(rna);

  for (i=0;i<n-1-THRESHOLD;i++) {/* leave room for j */
    if (secStr[i]=='.') 
      for (j=i+THRESHOLD+1;j<n;j++) {
	if (secStr[j]=='.') {
	  if (BP[i*n+j] && canAddBasePair(i,j,n,bps)) {
	    (*basePairsWhichCanBeAdded)[2*k]   = i;
	    (*basePairsWhichCanBeAdded)[2*k+1] = j;
	    k++;
	  } 
	}
      }
  }
  (*basePairsWhichCanBeAdded)[2*k] = -1;
}
/* Only base pairs that are present in the target structure can be added (analogous to GO model of folding for proteins)*/
void listOfNextStatesByAddingCompatibleBasePair(char *rna, char *secStr, int *bps, 
				      int **basePairsWhichCanBeAdded) {
  int i,j,k = 0;
  int n = strlen(rna);

  for (i=0;i<n-1-THRESHOLD;i++) {/* leave room for j */
    if (secStr[i]=='.') 
      for (j=i+THRESHOLD+1;j<n;j++) {
	if (secStr[j]=='.') {
	  if (BP[i*n+j] && canAddBasePair(i,j,n,bps) && targetStr[i]=='(' && targetStr[j]==')') {
	    (*basePairsWhichCanBeAdded)[2*k]   = i;
	    (*basePairsWhichCanBeAdded)[2*k+1] = j;
	    k++;
	  } 
	}
      }
  }
  (*basePairsWhichCanBeAdded)[2*k] = -1;
}

	/* 
	* Returns time increment to move to next state. Additionally,
	* update secStr and bps to next configuration by roulette wheel.
	*/
void nextConfig(char *rna, struct stateConfig * config, int *bpsa){
	/* 
	* Returns time increment to move to next state. Additionally,
	* update secStr and bps to next configuration by roulette wheel.
	*/
	float E0 = config->MFE;
	int ma = length(bpsa); /* Number of base pairs that can be added */
	int mr = length(config->bps);  /* Number of existing base pairs (can be removed) */
	char *tmp;    /* The temporary secondary structure */
	int i,j,k,l,x,n=strlen(rna);
	float E1,prob,flux;
	double z,timeInc,cumSum;
	
	int  *tmpbpsa,*tmpbps,tmpma,tmpmr;
	tmpbpsa = (int *) xcalloc(2*n*(n-1)/2+1,sizeof(int));
	tmpbps = (int *) xcalloc(n+1,sizeof(int));
	
	double *Pa = (double *) xcalloc(4*n*(n-1)/2+1,sizeof(float));
	double *Pr = (double *) xcalloc(4*n+1,sizeof(float));
	
	/*-------------------------------------------------------
	 Pa is prob of adding base pair (i,j), where
		Pa[4*k]   = 1/N(config->secStr) * min(1,exp(-(E(secStrBis)-E(config->secStr))/RT))
		Pa[4*k+1] = i
		Pa[4*k+2] = j
	 and to mark end of list of base pairs that can be added, we have
		Pa[4*k]   = -1
		Pa[4*k+1] = -1
		Pa[4*k+2] = -1  
	 Pr is prob of removing base pair (i,j), where format is identical to Pa. 
	-------------------------------------------------------*/
	tmp = (char *) xcalloc(n+1,sizeof(char)); //temporary secondary structure
	/* Copy the list of base pairs in bps into tmp. Note config->secStr[n]='\0'.*/
	for (i=0;i<=n;i++)    
		tmp[i] = config->secStr[i];
	flux = 0;
	for (x=0;x<ma;x++){ //compute tmp sec str by adding base pair (i,j) 
		i = bpsa[2*x]; j = bpsa[2*x+1];
		//k = length(config->bps);
		tmp[i] = '('; tmp[j] = ')';
		E1 = energy_of_struct(rna,tmp);
		if (DEBUG) printf("add base pair (%d,%u), MFE after adding:%lf\n",i,j,E1);

		if (!strcmp(method,"emc")){
			if (hasting==0){
				prob      = 1.0/(ma+mr) * MIN(1,exp(-(E1-E0)/RT));
			}
			else{
				k=0;
				while (config->bps[2*k]!=-1) {
					tmpbps[2*k] = config->bps[2*k];
					tmpbps[2*k+1] = config->bps[2*k+1];
					k++;
				}
				tmpbps[2*k] = config->bps[2*k];
				tmpbps[2*k+1] = config->bps[2*k+1];
				tmpbps[2*mr]=i; tmpbps[2*mr+1]=j; 
				tmpbps[2*(mr+1)]=-1; tmpbps[2*(mr+1)+1]=-1; //place new marker for end of base pairs
				listOfNextStatesByAddingCompatibleBasePair(rna,tmp,tmpbps, &tmpbpsa);
				tmpma = length(tmpbpsa);
				tmpmr = mr +1;
				prob = 1.0/(ma+mr) * MIN(1,exp(-((E1-E0)/RT)*((ma+mr)/(tmpma+tmpmr))));
			}
		}
		else if (!strcmp(method,"gil")){ //flux, not prob
			if (hasting==0){
				prob      = MIN(1,exp(-(E1-E0)/RT));
			}
			else{
				k=0;
				while (config->bps[2*k]!=-1) {
					tmpbps[2*k] = config->bps[2*k];
					tmpbps[2*k+1] = config->bps[2*k+1];
					k++;
				}
				tmpbps[2*k] = config->bps[2*k];
				tmpbps[2*k+1] = config->bps[2*k+1];
				tmpbps[2*mr]=i; tmpbps[2*mr+1]=j; 
				tmpbps[2*(mr+1)]=-1; tmpbps[2*(mr+1)+1]=-1; //place new marker for end of base pairs
				listOfNextStatesByAddingCompatibleBasePair(rna,tmp,tmpbps, &tmpbpsa);
				tmpma = length(tmpbpsa);
				tmpmr = mr +1;
				prob = MIN(1,exp(-((E1-E0)/RT)*((ma+mr)/(tmpma+tmpmr))));
			}
		}
		else{
			printf("error in method name\n");
			exit(1);
		}
		//if (DEBUG) printf("MFE:%lf,E1:%lf,RT:%lf,ma:%d,mr:%d,RT%lf\n",E0,E1,RT,ma,mr,RT);
		flux     += prob;
		Pa[4*x]   = prob; Pa[4*x+1] = i; Pa[4*x+2] = j; Pa[4*x+3] = E1;
		//now reset tmp to value before addition of (i,j)
		tmp[i] = '.'; tmp[j] = '.'; 
	}
	//indicate end of base pairs that can be added, not necessary since know ma
	Pa[4*x]   = -1; Pa[4*x+1] = -1; Pa[4*x+2] = -1; Pa[4*x+3]=-1;
	
	for (x=0;x<mr;x++){ //compute tmp sec str by removing base pair (i,j) 
		i = config->bps[2*x]; j = config->bps[2*x+1];
		tmp[i] = '.'; tmp[j] = '.';
		E1 = energy_of_struct(rna,tmp);
		if (DEBUG) printf("remove base pair (%d,%d), MFE after removing:%lf\n",i,j,E1);
		if (!strcmp(method,"emc")){
			if(hasting==0){
				prob      = 1.0/(ma+mr) * MIN(1,exp(-(E1-E0)/RT));
			}
			else{
				k = 0;
				while (config->bps[2*k]!=-1) {
					tmpbps[2*k] = config->bps[2*k];
					tmpbps[2*k+1] = config->bps[2*k+1];
					k++;
				}
				tmpbps[2*k] = config->bps[2*k];
				tmpbps[2*k+1] = config->bps[2*k+1];			
				/* Find position k of (i,j) in list, delete config->bps[2*k]=i and config->bps[2*k+1]=j
				by moving down all entries in bps. */
				k = 0;
				while (tmpbps[2*k]!=-1) {
					if (tmpbps[2*k]==i && tmpbps[2*k+1]==j)
					break;
					k++;
				}
				//move up entries in bps. Only do this when update secStr
				if (tmpbps[2*k]!=-1) {
					l = k + 1;
					while (tmpbps[2*l]!=-1) {
						tmpbps[2*(l-1)] = tmpbps[2*l];
						tmpbps[2*(l-1)+1] = tmpbps[2*l+1];
						l++;
					}
					tmpbps[2*(l-1)]=-1;
				}
				listOfNextStatesByAddingCompatibleBasePair(rna,tmp,tmpbps, &tmpbpsa);
				tmpma = length(tmpbpsa);
				tmpmr = mr - 1;
				prob = 1.0/(ma+mr) * MIN(1,exp(-((E1-E0)/RT)*((ma+mr)/(tmpma+tmpmr))));
			}
		}
		else if (!strcmp(method,"gil")){//flux used instead of prob
			if(hasting==0){
			prob      = MIN(1,exp(-(E1-E0)/RT));
			}
			else{
				k = 0;
				while (config->bps[2*k]!=-1) {
					tmpbps[2*k] = config->bps[2*k];
					tmpbps[2*k+1] = config->bps[2*k+1];
					k++;
				}
				tmpbps[2*k] = config->bps[2*k];
				tmpbps[2*k+1] = config->bps[2*k+1];			
				/* Find position k of (i,j) in list, delete config->bps[2*k]=i and config->bps[2*k+1]=j
				by moving down all entries in bps. */
				k = 0;
				while (tmpbps[2*k]!=-1) {
					if (tmpbps[2*k]==i && tmpbps[2*k+1]==j)
					break;
					k++;
				}
				//move up entries in bps. Only do this when update secStr
				if (tmpbps[2*k]!=-1) {
					l = k + 1;
					while (tmpbps[2*l]!=-1) {
						tmpbps[2*(l-1)] = tmpbps[2*l];
						tmpbps[2*(l-1)+1] = tmpbps[2*l+1];
						l++;
					}
					tmpbps[2*(l-1)]=-1;
				}
				listOfNextStatesByAddingCompatibleBasePair(rna,tmp,tmpbps, &tmpbpsa);
				tmpma = length(tmpbpsa);
				tmpmr = mr - 1;
				prob = MIN(1,exp(-((E1-E0)/RT)*((ma+mr)/(tmpma+tmpmr))));
			}
			
		}
		else{
			printf("error in method name\n");
			exit(1);
		}
		flux     += prob;
		Pr[4*x]   = prob; Pr[4*x+1] = i; Pr[4*x+2] = j; Pr[4*x+3] = E1;
		//now reset tmp to value before deletion of (i,j)
		tmp[i] = '('; tmp[j] = ')'; 
	}
	//Now normalize Pa,Pr by dividing by flux. This renders Pa,Pr cond prob
	for (x=0;x<ma;x++){
		Pa[4*x] = Pa[4*x]/flux;
		//printf("flux%d adding(%u,%u): %lf\n",x,Pa[4*x+1],Pa[4*x+2],Pa[4*x]);
		}
	for (x=0;x<mr;x++){
		Pr[4*x] = Pr[4*x]/flux;	
		//printf("flux%d removing(%u,%u): %lf\n",x,Pr[4*x+1],Pr[4*x+2],Pr[4*x]);
		}
	//Now use roulette wheel to determine next move
	z      = drand48();     //z in (0,1)
	//if (DEBUG) printf("z : %lf\n",z);
	cumSum = 0;
	for (x=0;x<ma+mr;x++){
		if (x<ma)
		cumSum += Pa[4*x];
	else
		cumSum += Pr[4*(x-ma)];
	//if (DEBUG) printf("cumsum : %lf\n",cumSum);
	if (z<cumSum) break;
	}
	
	timeInc = sampleExpDist(1/flux);
	if(verbose) printf("%s\t%lf\t%lf\t%d\n",config->secStr,config->MFE,timeInc,ma+mr);
	config->time += timeInc;
	//update the state configuration
	if (x<ma){
		i=Pa[4*x+1]; j=Pa[4*x+2];
		config->secStr[i]='('; config->secStr[j]=')';
		config->MFE = Pa[4*x+3];
		if (DEBUG) printf("next secStr : %s, MFE: %lf\n",config->secStr,config->MFE);
		config->bps[2*mr]=i; config->bps[2*mr+1]=j; 
		//now place new marker for end of base pairs
		mr += 1; config->bps[2*mr]=-1; config->bps[2*mr+1]=-1; 
	}
	else{ //ma<=x<ma+mr
		i=Pr[4*(x-ma)+1]; j=Pr[4*(x-ma)+2];
		config->secStr[i] = '.'; config->secStr[j] = '.';
		config->MFE = Pr[4*(x-ma)+3];
		if (DEBUG) printf("next secStr : %s, MFE: %lf\n",config->secStr,config->MFE);
		/* Find position k of (i,j) in list, delete config->bps[2*k]=i and config->bps[2*k+1]=j
		   by moving down all entries in bps. */
		k = 0;
		while (config->bps[2*k]!=-1) {
			if (config->bps[2*k]==i && config->bps[2*k+1]==j)
			break;
			k++;
		}
		//move up entries in bps. Only do this when update secStr
		if (config->bps[2*k]!=-1) {
			l = k + 1;
			while (config->bps[2*l]!=-1) {
				config->bps[2*(l-1)] = config->bps[2*l];
				config->bps[2*(l-1)+1] = config->bps[2*l+1];
				l++;
			}
			config->bps[2*(l-1)]=-1;
		}
		else
			printf("There is something wrong, bps too short?\n");
	}
	if (!strcmp(config->secStr,targetStr)) //update the number of times passing the target structure
		config->folded += 1 ;
	config->num      += 1; 
	if (DEBUG) {
		printSecStrListAsString(rna);
		printSecStrListAsString(config->secStr);
		printf("MFE:%lf\n",config->MFE);	
		printBasePairList(config->bps,sizeof(config->bps));
	}
	free(tmp); free(Pa); free(Pr);
	//inline function for  -1.0/flux*log(drand48());
}

struct stateConfig monteCarlo(char *rna) {
	int ma,mr,n = strlen(rna);
	struct stateConfig out;
	/* basePairsWhichCanBeAdded */
	int  *bpsa = (int *) xcalloc(2*n*(n-1)/2+1,sizeof(int));
	out.num         = 0; //total num of moves
	out.time        = 0; //total time for event driven MC procedure
	out.secStr      = initialEmptySecStr(rna); // initially empty 
	out.MFE		    = 0;
	out.bps         = (int *) xcalloc(n+1,sizeof(int));
	out.bps[0]      = -1;
	out.folded      = 0; //target str is not yet found
	
	if(!strcmp(out.secStr,targetStr)){
		printf("No trajectory found! target structure is the empty structure!\n");
		exit(0);
	}
	/* Update the list of base pairs that can be added, bpsa. */
	
	while (out.num < trjLength && trjLength!=-1 ) // when trajectory length is given get the path of length trjLength. good for getting long trajectories
	{
		if (GO)
			listOfNextStatesByAddingCompatibleBasePair(rna, out.secStr,out.bps, &bpsa);
		else
			listOfNextStatesByAddingBasePair(rna, out.secStr,out.bps, &bpsa);
		nextConfig(rna,&out,bpsa);	
	}
	while (strcmp(out.secStr,targetStr) && trjLength==-1 ) //stop at the target structure when trajectory length is not given
	{
		if (GO)
			listOfNextStatesByAddingCompatibleBasePair(rna, out.secStr,out.bps, &bpsa);
		else
			listOfNextStatesByAddingBasePair(rna, out.secStr,out.bps, &bpsa);
		nextConfig(rna,&out,bpsa);
	}
	
	if(!strcmp(out.secStr,targetStr)){
		if (GO)
			listOfNextStatesByAddingCompatibleBasePair(rna, out.secStr,out.bps, &bpsa);
		else
			listOfNextStatesByAddingBasePair(rna, out.secStr,out.bps, &bpsa);
		int ma = length(bpsa); 
		int mr = length(out.bps);
		if(verbose) printf("%s\t%lf\t%lf\t%d\n",out.secStr,out.MFE,0.0,ma+mr);
	}
	free(bpsa);
	return out;
}

int run(char *rna, int NUMRUNS) {
	struct timeval t;
	int n, i, folded;
	double firstPassageTimeList[NUMRUNS];
	double logfirstPassageTimeList[NUMRUNS];
	struct stateConfig mcout;
	double mean0,stdev0,max0,min0;
	
	/* Base pairing matrix */
	BP = canBasePair(rna);
	
	mfeSS = (char*) xcalloc(strlen(rna)+1,sizeof(char));
	minE = fold(rna,mfeSS);
	printf("#MFE: %g\t%s\n",minE,mfeSS);
	if (targetStr==NULL)
		targetStr = mfeSS;
	printf("#Target structure: %s\n",targetStr);
	/* use time in seconds to set seed */
	gettimeofday(&t, NULL);
	if (synchronized)
		srand48(seed);  //used for exact comparison between Gillespie and MC
	else
		srand48(t.tv_usec);

	n = strlen(rna);
	
	printf("#THRESHOLD = %u, TEMP = %g\n",THRESHOLD,temperature+273.15);
	
	for (i=0;i<NUMRUNS;i++) {
		printf("# trajectory %d\n",i+1);
		mcout = monteCarlo(rna);
		printf("Time:%g\n",mcout.time);
		firstPassageTimeList[i] = mcout.time;
		folded                  = mcout.folded; //if MFE str found
		if (verbose)
			printf("Number of MFE structure visits: %d\n",folded); 
		free(mcout.secStr);
		free(mcout.bps);
	}
	printf("Statistics after %d runs\n",NUMRUNS);
	mean0  = mean(firstPassageTimeList, NUMRUNS);
	stdev0 = std(firstPassageTimeList, NUMRUNS);
	max0   = max(firstPassageTimeList, NUMRUNS);
	min0   = min(firstPassageTimeList, NUMRUNS);
	
	if(MFPT){
	fprintf(timeFile,"MFE:%g\tTEMP:%g\tENG:Turner%s\n",minE,temperature,engModel);
		for (i=0;i<NUMRUNS;i++)
			fprintf(timeFile,"%d\t",firstPassageTimeList[i]);
		}
	printf("Mean: %f\tStDev: %f\tMax: %f\tMin: %f\n",mean0,stdev0,max0,min0);
	return 0;
}

int main(int argc, char *argv[]) {
	int i;
	char outFileName[121] ,*rna,*filePrefix="rna", fastaName[121]=""; 
	if (argc>1){
		for(i=1;i<argc;i++){
			if (argv[i][0] == '-'){
				switch (argv[i][1]){
					case 's':   //sequence
						if(i+1<argc && argv[i+1][0]!='-')
							rna = argv[i+1];
						else{
							printf("\nerror in the input sequence");
							printUsage(argv[0]);
						}break;
					case 'c':   //target structure
						if(i+1<argc && argv[i+1][0]!='-')
							targetStr = argv[i+1];
						else{
							printf("\nerror in the input target structure");
							printUsage(argv[0]);
						}break;
					case 'f':  //fasta file
						if (i+1<argc && argv[i+1][0]!='-'){
							strcpy(fastaName,argv[i+1]); 
						}
						else{
							printf("\nerror in the input file");
							printUsage(argv[0]);
						}break;
					case 'm': //method
						if (i+1<argc && (!strcmp("emc",argv[i+1])|| !strcmp("tmc",argv[i+1]) || !strcmp("gil",argv[i+1])))
							strcpy(method,argv[i+1]);
						else{
							printf("\nerror in method selection");
							printUsage(argv[0]);
						}break;
					case 'd':  //dangles
						if(argc>i+1 && argv[i+1][0]!='-')
							dangles = atoi(argv[i+1]);
						else{
							printf("\nerror in dangles");
							printUsage(argv[0]);
						}break;
					case 't':  //temperature
						if(argc>i+1 && argv[i+1][0]!='-'){
							temperature = atoi(argv[i+1]);
							RT= 0.00198717 * (273.15 + temperature);
						}
						else{
							printf("\nerror in temperature");
							printUsage(argv[0]);
						}break;
					case 'e':  //energy model
						if(argc>i+1 && argv[i+1][0]!='-')
							strcpy(engModel,argv[i+1]);
						else{
							printf("\nerror in energy model");
							printUsage(argv[0]);
						}break;
					case 'l':  // number of trajectory states
						if(argc>i+1 && argv[i+1][0]!='-')
							trjLength = atoi(argv[i+1]);
						else{
							printf("\nerror in trajectory length");
							printUsage(argv[0]);
						}break;
					case 'r': //number of times algorithm is run for the given rna
						if(argc>i+1 && argv[i+1][0]!='-')
							numRun = atoi(argv[i+1]);
						else{
							printf("\nerror in trajectory length");
							printUsage(argv[0]);
						}break;
					case 'y':  // give a NON random, fixed seed
						if(argc>i+1 && argv[i+1][0]!='-'){
							seed = atoi(argv[i+1]);
							synchronized = 1;
						}
						else{
							printf("\nerror in sync flag");
							printUsage(argv[0]);
						}break;
					case 'o': //output name
						if(argc>i+1 && argv[i+1][0]!='-')
							filePrefix = argv[i+1];
						else{
							printf("\nerror in output prefix");
							printUsage(argv[0]);
						}break;
					case 'g': //GO model
						GO = 1;
						break;
					case 'a'://Hasting's trick
						hasting=1;
						break;
					case 'v': //print trajectory
						verbose = 1;
						break;
					case 'h':
						printHelp(argv[0]);
						break;
					default : 
						printUsage(argv[0]);
						break;	
				}
			}
		}
	}
	else{
		printUsage(argv[0]);
		
	}
	if (strcmp(fastaName,""))
		rna = getInputRNA(fastaName);
	
	if (!validNucleotideString(rna)) {
		printf("Invalid nucleotide string\n");
		exit(1);
	}
	if(DEBUG) printf("input RNA is : %s\n",rna);	
	sprintf(outFileName,"%s_firstPassageTimes",filePrefix);
	if (MFPT){
		timeFile = fopen(outFileName,"w");
		fprintf(timeFile,"#seq:%s\tlength:%u\t",rna,(unsigned)strlen(rna));
	}
	if (targetStr!=NULL){
		if (strlen(rna) != strlen(targetStr)){
			printf("Sequence and the target structure must have the same length!\n");
			exit(1);
			}
		}
	//read the energy parameter file
	char * execPath = getExecPath(argv[0]);
	if (!strcmp(engModel,"99")) //use turner99 energy model
		strcat(execPath,turner99); 	
	else if (!strcmp(engModel,"04"))
		strcat(execPath , turner04);
	else if (!strcmp(engModel,"07"))
		strcat(execPath , andronescu07);
	else
		printUsage(argv[0]);
	read_parameter_file(execPath);
	
	run(rna, numRun);
	
	if (DEBUG) printf("native run done\n");
	
	if (MFPT) fclose(timeFile);
	free(mfeSS);
	return 0;
}
