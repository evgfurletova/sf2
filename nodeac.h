#ifndef NODEAC_H
#define NODEAC_H

#include "maindata.h"
#include <stdint.h>
#include <list>
#include "nodeov.h"

						//Classes of descriptors of nodes of AC-trie
//The vertexes of the trie are prefixes of words from motif
//The trie has prefix  and suffix links. The pair (s,t) is prefix link iff s is maximal prefix of t (s is prefix parent of t); 
//(s,t) is suffix link iff t is maximal suffix of s that is prefix of words from motif (t is suffix parent of s)
//The nodes of trie can be devided on internal vertexes and leaves

//Descriptor of common parameters of all vertexes of the trie

class NodeAC{
	public:
	uint8_t sign;			 //number of symbol in the alphabet which is assigned with prefix link leading  to the vertex
	bool leaf;				 //1, if it is a leaf; o- otherwise
	NodeAC* RParent;		 // suffix link 


	static int NumACNodes;	  // number of nodes in the trie 
    static NodeAC* ACRoot;	   //root of AC trie

	///////METHODS///////////////
	NodeAC(); //constructor

	//To insert new node to the AC trie
	//Input parameters:
	//node - prefix parent of new node; 
	//sign - symbol (number in alphabet) on prefix link leading from node to new node
	static void InsertNode(NodeAC* &node, int sign, int Number);
	static int InsertWord(int* word, int num);
    
    // Creates suffix links of the trie, marks nodes of the trie that are overlaps 
	static void CreateTrie();			
	
	//Prints word in the motif by depth-first traversal of the trie
	//Input parameters:
	//node - current considered node during traversal, word - word corresponding to the node
	static void  PrintWords(std::ostream*ff, NodeAC* node, std::string word);
};		


//Descriptor of internal node of the trie
class InternAC:
	public NodeAC
{
	public:
	int num;						  // number of the same node in the overlap graph
	std::list<NodeAC*> LChilds;       // List of prefix links 
	bool main;		   	              // if this vertex is overlap then main is 1, else 0;                      (1б)

	/////////METHODS//////////////////
	InternAC();
    
    // clear an instance of ptr like desctructor
    static void Clear(InternAC * ptr);
};


//Descriptor of terminal node (leaf) of the trie
class LeafAC:
	public NodeAC
{
	public:
	unsigned int leafnum;	     	  //number of the leaf													(4б)			

	/////////METHODS///////////////////
	LeafAC();
};

#endif