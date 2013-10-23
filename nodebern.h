#ifndef NODEBERN_H
#define NODEBERN_H

#include "nodeov.h"
#include "maindata.h"

////Descriptor of node w of OvGraf in case of Bernoulli model 

class NodeBern:
	public NodeOv				//class contaning common for all models parameters of descriptor  
{
public:

	double BackProb;			//Prob(Back(w))
	double* DeepProbs;          //Deep probabilities Prob(Back(H*(w,r)) ), r - right deep nodes

	double  WordProb;			//if the node w is a right deep node: 
	                            // Prob(HH~(w)), where for all H from HH*(w), w=rpred(H) 
   
	double* FirstTemp;			// If w is a right deep node then for all   
								//p = 1,..,NOccur contains probabilities Prob(RE(n, p, H(w))) 
								//calculated in a current step n; 
								//else the array is empty.

	double* ProbMark;			// for all p = 1,..,NOccur; k = n-|w|,..,n 
								// contains probabilities Prob(R(k, p, w)) for a current step n; 

	///METHODS
	NodeBern();				//constructor
	~NodeBern();			//destructor

	static void ProbCalc(void); //computation of Pvalue
};

//AUXILARY LISTS
extern double** BDProbs; //Let w be processed node, n- current stage; depth - depth of the path leading to w
						//x_0,..x_depth - overlap prefixes of w, x_depth = w
						//BDPobs[k][p] = Prob(D(n-m+|x_k|,p+1,x_k)), p= 0,...,p0-1, k = 0, ..,depth 

#endif
