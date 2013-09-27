#ifndef STATEDATA_H
#define STATEDATA_H

#include <vector>

typedef struct PriorList
{
	int pos;				 // pointer to a state;
	double prob;			 // probability
}PriorList;


/////////Descriptor of pair <w,q>, w in OV(HH), q in AllState(w) (common part for HMM and Markov models)///////////////////////////
typedef struct H_M_State
{
	int ID;					// identificator of state q
	int RParentPos;			// number of the same state in list rpred(w).States

	int NumBackProbs;		//|PriorList(w,q)| (the notation is given in our paper about HMM)
	PriorList* BackProbs;   //  Back probabilities. For all q' in PriorList(w,q),  BackProbs[q'].pos is position of q' in lpred(w).States
							//BackProbs[q'].prob = Prob(q', Back, q)

	int* NumDeepProbs;		//list of |ReachState(q,H*(w,r))|, r is a right deep node such that exist class  H*(w,r)  (see paper)
    PriorList** DeepProbs;  //Deep probabilities. For all s in ReachState(q,H*(w,r)), DeepProbs[r][s] .pos is position of s in r.States
							//DeepProbs[r][s].prob = Prob(q, Back(H*(w,r)),s)

	double* FirstTemp;		//array to store temporary value Prob(RE(n,p,w),q), p =1,..,p0 
	double* ProbMark;		// for all p = 1,..,p0; k = n-|w|,..,n 
							// For all s in ReachState(q,H*(w,r)), contains probabilities Prob(R(k, p, w),q) for a current step n; 

							//WordProbs[q'].prob = Prob(q,H~(w),state)	

	 
	/////methods////////////
	H_M_State();
	void Clear(int NDLinks,int rflag); //clearing of data in descriptor
}H_M_State;




/////////Descriptor of pair <w,q>, w in OV(HH), q in AllState(w) (additional parameters for HMM )///////////////////////////
typedef struct HMM_State:
	public H_M_State
{
	int NumWordProbs;				//|ReachState(H~(w),q)| (see paper),H~(w) is thw set of words H from HH, such that r = rpred(H) 
	PriorList* WordProbs;			// Word probabilities. For all q' in ReachState(H~(w),q), WordProbs[q'].pos = q'
									//WordProbs[q'].prob = Prob(q,H~(w),state)	
		
	HMM_State();
}HMM_State;


#endif
