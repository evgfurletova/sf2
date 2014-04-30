#include "nodeac.h"
#include "nodeov.h"
#include "nodebern.h"
#include "h_m_node.h"
#include "statedata.h"
#include "m_trtree.h"
/////////////////NODES OF AC TRIE////////////////////////////////
                 /////////NodeAC//////////
int NodeAC::NumACNodes = 0;
NodeAC* NodeAC::ACRoot = nullptr;

NodeAC::NodeAC(){
//	ID = NodeAC::NumACNodes;
	RParent = nullptr;
}



                   /////InternAC/////

InternAC::InternAC()
	:NodeAC()
{
	leaf = 0;
	main = 0;
	num = 0;

}

void InternAC::Clear(InternAC * ptr){
	ptr->LChilds.clear();
};

					///LeafAC///////
LeafAC::LeafAC()
	:NodeAC()
{
	leaf = 1;
	leafnum = 0;
}

/////////////////NODES OF OVERLAP GRAPH////////////////////////////////

NodeOv* NodeOv::Root;			
    int NodeOv::NumOVNodes = 0;		    
    int NodeOv::NumRDNodes = 0;			
    int NodeOv::NumLDNodes = 0;		
    int NodeOv::NClasses = 0;	

NodeOv::NodeOv() {
	num = 0;
    len = 0;
	rdeep = 0;			  
	ldeep = 0;		    
    rootchild = false;
	RParent = nullptr;		
	NLChilds = 0;
	NRChilds = 0;	      	    
	NDLinks = 0;		
    LChilds = nullptr; 
	RChilds = nullptr;    
	DeepLinks = nullptr; 	
};		

NodeOv::~NodeOv(){		
    free(LChilds);   
	free(RChilds); 
	free(DeepLinks); 	
};		

NodeOv* NodeOv::NewNode(){
	NodeOv* newnode;
	if(MainData::order == 0){
		newnode = new NodeBern();
	}
	else{
		newnode = new H_M_Node();
	}

	return newnode;
};
/////////////////NODES OF OVERLAP GRAPH FOR BERNOULLI////////////////////////////////

NodeBern::NodeBern()
    : NodeOv()
{
    BackProb = 0;
    DeepProbs = nullptr;
    WordProb = 0;
    rootchild = false;
    FirstTemp = nullptr;
    ProbMark = nullptr;
};

NodeBern::~NodeBern(){	
	free(DeepProbs);       
	free(FirstTemp);			
	free(ProbMark);		
};

/////////////////NODES OF OVERLAP GRAPH FOR HMM////////////////////////////////


H_M_Node::H_M_Node()
{
	NodeOv();
	States = nullptr;
	NStates = 0;
}


H_M_Node::~H_M_Node()
{
	int i; 
	for(i = 0; i < NStates; i++){
		States[i]->Clear(NDLinks,this->rdeep);
		delete States[i];
		States[i] = nullptr;
	}
	free(States);
	States = nullptr;

}


H_M_State::H_M_State()
{
	ID = -1;					
	RParentPos = -1;			

	NumBackProbs = 0;		
	BackProbs = nullptr;  
							

	NumDeepProbs = nullptr;
	DeepProbs = nullptr;
    FirstTemp = nullptr;		
	ProbMark = nullptr;					
}



HMM_State::HMM_State()
	:H_M_State()
{
	NumWordProbs = 0;				
    WordProbs = nullptr;			
}


void H_M_State::Clear(int NDLinks, int rflag)
{	
	free(BackProbs);   

	free(NumDeepProbs);	

    int i;
    for(i = 0; i < NDLinks; i++){
        delete[] DeepProbs[i];
		DeepProbs[i] = nullptr;
	}
    delete[] DeepProbs;
	DeepProbs = nullptr;
	free(FirstTemp);
	free(ProbMark);	
	FirstTemp = nullptr;
	ProbMark= nullptr;

	if((MainData::order <0)&&(rflag == 1)){
		HMM_State* state = static_cast<HMM_State*>(this);
		free(state->WordProbs);				
	}

}

/////////Nodes of MTrTRee/////////////
M_TrTree::M_TrTree(){
	sign = -1;	
	len = 0;			
	Childs = nullptr;
	NStates = 0;
	States = nullptr;
};


M_TrTree::~M_TrTree(){
	
	free(States);
    int i;
	for(i = 0; i < MainData::AlpSize; i++){
		if(Childs[i] != nullptr){
			delete Childs[i];
			Childs[i] = nullptr;
		}
	}

    delete[] Childs;
	Childs = nullptr;
 
};

KM_TrTree::KM_TrTree()
    : M_TrTree()
{
//	M_TrTree();
	NumWLinks = 0;
	WLinks = nullptr;
	NumWProbs = nullptr;
	WProbs = nullptr;
};
	
KM_TrTree::~KM_TrTree(){
	free(WLinks);
	free(NumWProbs);
	int i;
	for(i = 0; i < NumWLinks; i++){
        delete[] WProbs[i];
		WProbs[i] = nullptr;
	}
    delete[] WProbs;
	WProbs = nullptr;
};
	
