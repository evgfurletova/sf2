//#include "StdAfx.h"
#include "nodebern.h"

							//NOTATIONS
		//Len p0 - be desired number of occurences; n0- text length; 
		//n - current stage of the algorithm; m - pattern length; Alp - alphabet; HH - pattern

double*	BSumProb = nullptr;	     //BSumProb[j] = Sum_{k<=n}R(k,j,HH)

double** BDProbs = nullptr;		//Let w be processed node, n- current stage; depth - depth of the path leading to w
						//x_0,..x_depth - overlap prefixes of w, x_depth = w
						//BDPobs[k][p] = Prob(D(n-m+|x_k|,p+1,x_k)), p= 0,...,p0-1, k = 0, ..,depth 

//////////////////////////////////////
// Initialization of lists ProbMark and FirstTemp in descriptors of nodes of overlap graph

void IniProbs(NodeBern* node){
		int i, row;
	
		row = MainData::WordLen - node->len;  

		node->ProbMark = Malloc<double>(MainData::NOccur*row);
		node->FirstTemp = Malloc<double>(MainData::NOccur);

		for( i = 0; i < node->NLChilds; i++){
			NodeBern* LC = static_cast<NodeBern*>(node->LChilds[i]);
			IniProbs(LC);
		}
		return;
	}


//calculate initial probabilities Prob(R(m,p,H(w))) and add it to w.ProbMark, 
//by bottom-up traversal of overlap graph following by right links

//Input parameters: node - processing node w
void Stepm(NodeBern* node){
	int i;
    int pos = MainData::WordLen - node->len - 1;

    if(node->rdeep == 1){								//if node is right deep
		node->ProbMark[pos] = node->WordProb;			//Add Prob(H~(w))(subset of words from HH, having right predecessor w) 
														//to ProbMark (stored in parameter WordProb)
	}
	for ( i = 0; i < node->NRChilds; i++){
		NodeBern* RC = static_cast<NodeBern*>(node->RChilds[i]);
		Stepm(RC);										//Run the procedure for right successors of  w
        int pos1 = MainData::WordLen - RC->len - 1;	
		node->ProbMark[pos]  += RC->ProbMark[pos1];		//Prob(R(m,p,H(w))) = Prob(RE(m,p,H(w))) + SumProb_{t,rpred(t) = w}(R(m,p,H(t))) 
	}
	
	if((node == NodeOv::Root)&&(MainData::NOccur == 1)){
		MainData::Pvalue = node->ProbMark[pos];
	}
}


//Computing of probabilities Prob(RE(n,p,H(r))), r is right deep node
//Input parameters: node - processing node w;
//n - Current stage minus length m of words in the motif (will think that n is stage of algorithm)
//depth - number of nodes on the path leading from root to w
void RProbCalc(NodeBern* node, int n, int depth){
	int i,j;

	int row = MainData::WordLen - node->len;
	int pos = n % row;							//position in r.ProbMark, contaning  Prob(R(n-m+|w|,p,H(w)))
												//at the and of the stage n it contains Prob(R(n,p,H(w)))
	memset(BDProbs[depth], 0x00, MainData::NOccur * sizeof(double)); //initialization of BDProbs[depth]
	
	//Internal loop. j = 0,...,p0 -1 
	for(j = 0; j < MainData::NOccur; j++){
		if(node != NodeOv::Root){
			/////////1. Computation of D-sets probabilities Prob(D(n-m+|w|,j,r))//////////////////////////
			if(node->rootchild == 1){ //if lpred(w) == root
				BDProbs[depth][j] = node->ProbMark[j*row + pos];									    // Prob(D(n-m+|w|,j,w)) = Prob(R(n-m+|w|,j,w))
			}
			else{
				BDProbs[depth][j] = BDProbs[depth -1][j]*node->BackProb + node->ProbMark[j*row + pos];  // Prob(D(n-m+|w|,j,w)) = Prob(D(n-m+|w|,j,w))*Prob(Back(w))+ Prob(R(n-m+|w|,j,w))
			}

		
		////////2. Computation of RE-sets probabilities////////////////// ////////	
	
		//////////A. Computation of Prob(C(n,p-1,r))-Prob(C(n,p,r)) for all r such that exist a class h*(w,r)
			if(node->ldeep == 1){				//if w is left deep node
				for(i = 0; i < node->NDLinks; i++){
				    NodeBern* rnode = static_cast<NodeBern*>(node->DeepLinks[i]);
					if(j == 0){ 
						rnode->FirstTemp[0] -= BDProbs[depth][0]* node->DeepProbs[i];	 // -Prob(D(n-m+|w|,j,w))*Back(h*(w,r))  
					}
					else{
						rnode->FirstTemp[j] += (BDProbs[depth][j-1] - BDProbs[depth][j])* node->DeepProbs[i];  //+(Prob(D(n-m+|w|,j-1,w))-Prob(D(n-m+|w|,j,w)))*Back(h*(w,r))  
					}
				}
			}
		}
	
		//////////A. Computation of Prob(F(n,p-1,w))-Prob(F(n,p,w)), if w is right deep node
		if(node->rdeep == 1){
			if(j == 0){
				node->FirstTemp[0] += (1 - BSumProb[0])* node->WordProb;	//+(1-Sum_{k<=n-m} R(k,1,HH))*Prob(H~(w))
			}
			else{
				node->FirstTemp[j] += (BSumProb[j-1] - BSumProb[j])*node->WordProb; //+Sum_{k<=n-m} (R(k,j-1,HH) - R(k,j,HH))*Prob(H~(w))
			}	
		}
	}
	////////////////////////Recursion, processing of left successors of w//////////////////////////////////////////////
	if(node != NodeOv::Root) 
		depth++;

	for ( i = 0; i < node->NLChilds; i++){
		NodeBern* LC = static_cast<NodeBern*>(node->LChilds[i]);
		RProbCalc(LC, n,depth);
	}
	return;
}



//Updating ProbMark by bottom up traversal of  OVGraf following suffix links. Computing  Prob(R(n,p,H(w)))
//Input parameters: node - current processing node of  OVGraf; 
void ROGUpdate(NodeBern* node, int n){
	int i,j;
	int row = MainData::WordLen - node->len;
	int pos = n % row;						//position in w.ProbMark to record Prob(R(n,p,H(w)))

	
	for(i = 0; i < node->NRChilds; i++){
		NodeBern* RC = static_cast<NodeBern*>(node->RChilds[i]);
		ROGUpdate(RC,n);							//processing of right successors of w

		for(j = 0; j < MainData::NOccur; j++){
			node->FirstTemp[j] += RC->FirstTemp[j]; //Prob(R(n,p,H(w))) = Prob(RE(n,p,H(w))) + Sum_{t,w=rpred(t)}Prob(R(n,p,H(t)))
			RC->FirstTemp[j] = 0;
		}	

	}
	
	for(j = 0; j < MainData::NOccur; j++){
		node->ProbMark[j*row + pos] = node->FirstTemp[j];  //record Prob(R(n,p,H(w))) to w.ProbMark
		if(node == NodeOv::Root){
			node->FirstTemp[j] = 0;
		}
	}



	return;
}


/////////////////DEBUG!!!////////////////////

void DebBernPrint(NodeBern* node, int step){
	int i,j;	

	if(node->num == 0){
		cout<<"Step "<<step<<"\n-------\n--------\n\n";
	}
	
	int row = MainData::WordLen - node->len;

	cout<<"Num of node: "<<node->num<<'\n';
	cout<<"ProbMark \n";
	for(i = 0; i< row; i++){ 
		for(j = 0; j < MainData::NOccur; j++){
			cout<<std::setprecision(15)<<node->ProbMark[j*row + i]<<'\t';
		}
			cout<<'\n';
	}
	cout<<"\n\n";
	for ( i = 0; i < node->NLChilds; i++){
		NodeBern* G = static_cast<NodeBern*>(node->LChilds[i]);
		DebBernPrint(G,step);
	}
}

/////////////////////////////////////


//Main function for p-value calculation

void NodeBern::ProbCalc(void){
	int n,j;
	///////1. Initialization//////////
	NodeBern* root = static_cast<NodeBern*>(NodeOv::Root);
	BSumProb  = Malloc<double>(MainData::NOccur);
	BDProbs = Malloc<double*>(MaxDepth); 
	for(j = 0; j < MaxDepth; j++)
		BDProbs[j] =  Malloc<double>(MainData::NOccur);

	IniProbs(root);
	
	///////2. Computation of initial probabilities of R-sets///////////////////
	if(MainData::TLen > MainData::WordLen-1){
		Stepm(root);
	}
	
	if(MainData::NOccur == 1){
		MainData::Pvalue = root->ProbMark[MainData::NOccur*MainData::WordLen - 1];
	}

	///////3. Main loop////////////

//	DebBernPrint(root,MainData::WordLen);
	if(MainData::TLen > MainData::WordLen){
		int i = 0;
		for(n = MainData::WordLen+1; n<= MainData::TLen; n++){	
			
			///n -th stage ////////
			int pos = i % MainData::WordLen;

			
			for(j = 0; j < MainData::NOccur; j++){	
				BSumProb[j] = BSumProb[j] + root->ProbMark[j* MainData::WordLen + pos]; // Compute probabilities of B-sets
			}
		
			RProbCalc(root, i,0);		// Computing probabilities of RE-sets
		    ROGUpdate(root, i);			//Updating data in graph. Computing probabilities of R-sets
			

			
			MainData::ProbRes = root->ProbMark[(MainData::NOccur - 1)*MainData::WordLen + pos];
			MainData::Pvalue += MainData::ProbRes; //Pvalue
			i++;
		}
	}
	free(BSumProb);
	for(j = 0; j < MaxDepth; j++)
		free(BDProbs[j]);
	free(BDProbs);

	NodeOv::DeleteGraf(NodeOv::Root);
	return;
}



