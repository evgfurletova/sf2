#include "nodeac.h"
#include "nodeov.h"



//To insert new node to the AC trie
//Input parameters:
//node - prefix parent of new node; 
//sign - symbol (number in alphabet) on prefix link leading from node to new node

void NodeAC::InsertNode(NodeAC* &node, int sign, int Number)
{	
	NodeAC* newnode;
	InternAC*  node1 = static_cast<InternAC*>(node);
	
	//1. Creating of new node
	NodeAC::NumACNodes ++;
	if (Number == 0) //if node is internal
	{
		newnode = new InternAC();	
	}
	else //if node is leaf 
	{ 
		newnode = new LeafAC();
		LeafAC*  newleaf = static_cast<LeafAC*>(newnode);
		newleaf->leafnum = Number;	
	}

	//2. Creating of prefix link leading to new node 
	newnode->sign = sign;	
	node1->LChilds.push_back(newnode);
	node = newnode;
	return;
}

int NodeAC::InsertWord(int* word, int num){
	int res = 0;
	InternAC* node =static_cast<InternAC*>(NodeAC::ACRoot); 
	int i;
	std::list<NodeAC*>::iterator j;
	for(i = 0; i < MainData::WordLen; i++){
		j = node->LChilds.begin();
		int flag = 0;
		while((flag == 0) &&(j != node->LChilds.end())){
			NodeAC* node1 = *j;
			if(node1->sign == word[i]){
				flag = 1;
				node = static_cast<InternAC*>(node1);
			}
			j++;
		}
		if(flag == 0){
		    NodeAC* node1 = node;
			if(i< MainData::WordLen - 1){
				InsertNode(node1,word[i], 0);
				node = static_cast<InternAC*>(node1);
			}
			else{
				InsertNode(node1,word[i], num);
				res = 1;
			}
		}
	}
	return res;
}

void preLink()
{
	NodeAC** queue = Malloc<NodeAC*>(NodeAC::NumACNodes + 1); //deque

	int first = 0;		//first element of deque
	int last = -1;      // last element of deque
	int size = 0;		//size of deque

	
	/////////1. Creating suffix links for nodes of length 1(prefix childs of root) ///
	
	//////// Add childs of root into the queue ///////
	InternAC*  root = static_cast<InternAC*>(NodeAC::ACRoot); // root of the trie

	root->main = 1;
	
	std::list<NodeAC*>::iterator i,j;
	
	for(i = root->LChilds.begin(); i != root->LChilds.end(); i++){
			NodeAC* LCnode = *i;
		    LCnode->RParent = root; //create suffix links for childs of root 

			last ++;
			size ++;
			queue[last] = LCnode; // add childs of root to the deque 
	}	

	while (size > 0) //queue is not empty
	{
	
	//////////2. Processing nodes from queue//////////////////////
		NodeAC *v1= queue[first]; //pop	first element of deque			
		
		InternAC* node = static_cast<InternAC*>(v1);
		first++;
		size--;
	
		for(i = node->LChilds.begin(); i != node->LChilds.end(); i++){
			//2.1. Consider prefix child LCnode of current node from deque
		    	NodeAC *LCnode = *i;	 //next prefix child of node LCnode
				int sign = LCnode->sign; //sign of edge leading from node to LCnode
				
				if(LCnode->leaf == 0){ 
					last ++;
				    size ++;
					queue[last] = LCnode; //push child of node to the deque
				}
				

				InternAC* w = static_cast<InternAC*>(node->RParent); //maximal suffix-overlap of node
				
			//2.2. finding of maximal suffix-overlap of LCnode
				
				int flag = 0;
				NodeAC* RPnode; 
				

				while (flag == 0){
					j = w->LChilds.begin();
					while((flag == 0)&&(j != w->LChilds.end())){  
					   NodeAC* v = *j;
					   if(v->sign == sign){
						   flag = 1;
						   RPnode = v;
					   }
					  j++;
					}

					if(flag == 0){
						if(w == root){
							flag = 1;
							RPnode = root;
						}
						else w = static_cast<InternAC*>(w->RParent);
					}

				}	/*while (flag == 0)*/

				LCnode->RParent = RPnode;
				
			}/*for(i = node->LChilds.begin(); i != root->LChilds.end(); i++)*/
	}/*while (size > 0)*/
	free(queue);
	return;
}

////////////////////////////////////////////////

///////////////////////////////////////

// It markes overlap nodes by traversal of suffix links, starting with a leaf
//Input Parameters: node - current node
void MainNode(NodeAC* node){
	while(node != NodeAC::ACRoot){
		node = node->RParent;
		
		InternAC* node1 = static_cast<InternAC*>(node);
		
		if(node1->main == 1){
			return;
		}
		node1->main = 1;
	}
	return; 
} 

//It markes all overlap nodes
void  MarkMainNodes(NodeAC* node){
	
	if(node->leaf == 1){
		MainNode(node);
	}
	else{
		InternAC* node1 = static_cast<InternAC*>(node);
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			MarkMainNodes(node2);
		}
	}
}

//Assigning of numbers to overlap nodes by depth-first traversal of the Ac trie
void  NumMainNodes(NodeAC* node){
	
	if(node->leaf == 0){
		InternAC* node1 = static_cast<InternAC*>(node);
		if(node1->main == 1){
			node1->num = NodeOv::NumOVNodes;
			NodeOv::NumOVNodes ++;
		}
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			NumMainNodes(node2);
		}
	}
}

int nw = 1;
void ReverseOrder(NodeAC* node){
	if(node->leaf == 0){
		std::list<NodeAC*>::iterator i;
		int k;

		InternAC* node1 = static_cast<InternAC*>(node);
		NodeAC* node2;
		vector<NodeAC*> Childs;
		Childs.resize(MainData::AlpSize);
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			Childs[node2->sign] = node2;
		}
		i = node1->LChilds.begin();
		for(k = 0; k < MainData::AlpSize; k++){
			if(Childs[k] != nullptr){
				*i = Childs[k];
				i++;
			}
		}

		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			ReverseOrder(*i);
		}
	}
	else{
		LeafAC* node1 = static_cast<LeafAC*>(node);
		node1->leafnum = nw;
		nw++;
	}
}


//////////////////////////////////////////

	//FOR DEBUG!!!//

void  DebTriePrint(NodeAC* node, string word){
	string word1 = word += MainData::AlpMas[node->sign];
	if(node->leaf ==0){
		InternAC* node1 = static_cast<InternAC*>(node);
		if(node1->main == 1){
			cout<<" Num: "<<node1->num<<'\t'<<word1<<'\n';
		}
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			DebTriePrint(node2,word1);
		}
	}
}

void  DebTriePrint1(NodeAC* node, string word){
	string word1 = word + MainData::AlpMas[node->sign];
	InternAC* RP = static_cast<InternAC*>(node->RParent);
//	cout<<"Node: "<<word1<<" Num: "<<node->ID;
	//cout<<'\n'<<"LPred: "<<word<<" RPred: "<<RP->ID<<"\n";
	
	if(node->leaf ==0){
		
		InternAC* node1 = static_cast<InternAC*>(node);
		cout<<"Main: "<<node1->main<<" OvNum: "<<node1->num<<"\n\n\n";
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			DebTriePrint1(node2,word1);
		}
	}
	else{
		cout<<"\n\n\n";
	}
}


void  DebTriePrint2(NodeAC* node, string word){
	string word1 = word += MainData::AlpMas[node->sign];
	if(node->leaf ==0){
		InternAC* node1 = static_cast<InternAC*>(node);
		cout<<" Num: "<<node1->num<<'\t'<<word1<<'\n';
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			DebTriePrint2(node2,word1);
		}
	}
}

/////////////////////////////////////////////////////
//Prints word in the motif by depth-first traversal of the trie
//Input parameters:
//node - current considered node during traversal, word - word corresponding to the node

void  NodeAC::PrintWords(ostream*ff, NodeAC* node, string word){
	if(node != NodeAC::ACRoot){
	     word.push_back(MainData::AlpMas[node->sign]);
	}
	if(node->leaf != 0){
		LeafAC* leaf = static_cast<LeafAC*>(node);
		*(ff)<<leaf->leafnum<<".  "<<word<<'\n';
	}
	if(node->leaf ==0){
		InternAC* node1 = static_cast<InternAC*>(node);
		std::list<NodeAC*>::iterator i;
		for(i = node1->LChilds.begin(); i != node1->LChilds.end(); i++){
			NodeAC* node2 = *i;
			PrintWords(ff,node2,word);
		}
	}
}

//////////////////////////////////////////
 // Creates suffix links of the trie, marks nodes of the trie that are overlaps 
void NodeAC::CreateTrie(void){
	if(MainData::mode < 2) ReverseOrder(NodeAC::ACRoot);
	//NodeAC::PrintWords(&std::cout, NodeAC::ACRoot,"");
	preLink();					     //creating of suffix links
	MarkMainNodes(NodeAC::ACRoot);   //marking of overlap nodes
	NumMainNodes(NodeAC::ACRoot);	 //numbering of overlap nodes
	//DebTriePrint2(NodeAC::ACRoot,"");
}
////////////////////////
