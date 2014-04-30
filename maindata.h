#ifndef MAINDATA_H
#define MAINDATA_H

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <deque>
#include "type_malloc.h"
#include "nodeac.h"
#include "sf_main.h"
//#include "nodeov.h"
using namespace std;

#define		FileLineMax		1000000
#define		AlpMax			100
#define		WordLenMax		2000
#define		MaxOccur		1000
#define		NMaxStates		2000
#define		MaxSeqLen		300

class NodeAC;

class MainData
{
public:
	// 1. Input files
	static std::string AlpFileName;
	static std::string PatternFileName;
	static std::string PSSMFileName;
	static std::string FootPrintFileName;
	static std::string ConsAlpFileName;
	static std::string OutName;
	static std::string OutWordsName;
	static std::string SeqFileName;

	// 2. Output files
	static ofstream ResLOG;						// output files
	static ofstream ResROG;
	static ofstream ResF;
	static ofstream ResWords;
	

	// 3. Input Parameters
	// 3.1. Alphabet and model of letters distribution 
	static int			 AlpSize;					// size of the alphabet
	static char		     AlpMas[AlpMax];			//alphabet
	static int		 	 order;						//Order of Markovian model (0 in the case of Bernoulli model)
	static double		 BernProb[AlpMax];			//Array of Bernulli probabilities for letters

	// 3.2 Parameters for the pattern
	//3.2.1 common parameters for all type of input parameters respect to the pattern
	static int	mode;						//type of patter description
	static int NWords;						// number of words in pattern
	static int WordLen;						// length of the pattern
	static int CrDistribFlag;
	static int StrandType;					//1 - if complement motif is added, 0 -otherwise
	static int DistribType;					//1 - if model parameters have to be estimated, 0 - if the parameters are given
	//3.2.2. parameters for the pattern described by PSSM

	static double**	PssmMas;				// Pssm matrix
	static double** RPssmMas;
	static double Thr;						// Threhold for Pssm
	//3.2.3 parameters for the pattern described by a motif, number of replacements and number of constant positions
	static int* motif;						//motif
	static std::string motifstr;			//motif
	static int Nreplace;					//number of replacements in the pattern
	static int*	ConstPositions;				//list of constant positions		
	static int  NumConstPositions;
	//3.2.4 parameters of the pattern described by a consensus
	static int* consensus;					//consensus
	static std::string consensusstr;
	static vector< vector<char> > ConsAlp;	//alphabet for consensus
	//3.2.5 parameters for random pattern
    static double* RandPatProbs;			//probabilities to generate a random pattern
	// 3.3 Information for calculating of probabilities
	static int TLen;						//length of random text 
	static int NOccur;						// Maximal number of patterns occurrences in the text
	
	// 4. Variables
	static int MAX;
	// 4.1. Variables for graphs creating
	
	// 4.2	Variables for probabilities computation
	static double	ProbRes;				//Sum_q(Prob_q(R(TLen,NOccur,HH))) (see paper)
	static double	Pvalue;					//Probability to fined desired number of occurences in a random text 
	// 4.2.1. for a Bernoulli model

	//Let Q be a set of states 
	// 4.2.2. for a Markov model
	static double**	MarkovProbs;			//Matrix probabilities of symbols from the alphabet for different states;
											// For Markov model the states are in alphabetical order;
	static int MarkovType;					//1, if initial distribution is given; 0 - if Calculated 
	// 4.2.3. for a non-deterministic HHM
	static double*** ND_HHMProbs;			//Mairix (Q*Q*Alp->П(q,a,q')) of probabilities  П(q,a,q') being in a state q to generate a symbol a and traverse to a state q’.
	static vector<int>** ND_HHMTrans;		//Matrix (Q*Alp->Q)of sets of states such that П(q,a,q')>0
	// 4.2.4. for a deterministic HHM
	static double** D_HHMProbs;				//Mairix (Q*Alp->P(q,a)) of probabilities  P(q,a) being in a state q to generate a symbol a
	static int** D_HHMTrans;				//Mairix (Q*Alp->Ф(q,a)) state q=Ф(q,a) in which system passes from a state q generating a symbol a
	
	//5.
	static int Error;						//the code of error
	static int Out_mode;					//For the Markov model. if 1, then initial probabilities is given, else to be generated

	//6. AC trie information
	

	//7. Overlap graph


	//////////////////////////////



public:

	static	int		AToi(char Let);						//gives number of a symbol Let in AlpMas
	static	char    IToa(int pos);						//gives symbol Let having number pos in AlpMas
	static  int		in_Alp(char *str);					//checks is str a word in the alphabet
	static  int     EstimateDistribParams(void);		//estimates parameters of probabilities model
	static  int	    GenRanWords(void);					//generates random pattern
	static  double  CountThr(char *word);				//calculates cut-of for the pattern described by PSSM
	static	void    SetScorMas(double* SMas, double** Matrix);
	static	void	ReverseMatrix(void);
	static  int     GenPssmWords1(NodeAC* node, int *word, string& stword, int i, double score, double* SMas); //generates pattern described by PSSM
	static  int     GenPssmWordsR(NodeAC* node, int *word, string& stword, int i, double score, double scoreR, double* SMas, double* SMasR);
	static  void    MotifVariations1(NodeAC* node, int num, int pos, int* word);			//generates pattern described by a word and number of mismatches
	static  int     Pos_In_Cons_Alp(char Let);									//gives position of a letter Let in consensus alphabet
	static  int     ConsVariations1(NodeAC* node, int pos, int* word);			//generate pattern described by consensus
	static  int		GetInput(void);												//sets input parameters
	static	int		ComLineParse(int argc, char **argv);						//parsing of the command line
	static	void	CrOutFiles(void);											//prints information to the out files
	static	void    ErrorDetect(int Error);										//detects type of arised error
	static  void	PrintMain(ostream *ff,int exitflag);						//prints main report about results of program work
	static  int		Check_Out_of_Range(void);									//checks input parameters
};

//MainData Data;
extern std::string Words;
#endif
