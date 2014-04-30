#include "nd_hhm_prob.h"
#include "h_m_node.h"
#include "maindata.h"
#include <iostream>
#include <sstream>

double ND_HHM_Prob::TermCondProb(int q1,const std::string & word, size_t wordsize, size_t from, int q2){
    if(wordsize-from == 0){
		if(q1 == q2){
			return 1;
		}
		else{
			return 0;
		}
	}
	
	double prob = 0;
	size_t i;
    int pos = MainData::AToi(word.at(from));
	if(wordsize == 1){
		prob = MainData::ND_HHMProbs[q1][q2][pos];
	}
	else{
		size_t s = MainData::ND_HHMTrans[q1][pos].size();
//		word.erase(0,1);
		for(i = 0; i < s; i++){
			int state = MainData::ND_HHMTrans[q1][pos][i];
			double p = MainData::ND_HHMProbs[q1][state][pos];
            prob += TermCondProb(state,word,wordsize,from+1,q2)*p;
		}
	}
	return prob;
};


void ND_HHM_Prob::AllTermProbs(int q0, double* PredProbs, int* PredStates, int npredstates, double* NewProbs, int* NewStates, const std::string & word, size_t wordsize, size_t from){
	
	size_t i, s;
    int pos = MainData::AToi(word.at(from));
	int k,q1,q;
	double p;

	int nstates;
	if(from == 0){
		s = MainData::ND_HHMTrans[q0][pos].size();
		nstates = s;
		for(i = 0; i < s; i++){
			q = MainData::ND_HHMTrans[q0][pos][i];
			p = MainData::ND_HHMProbs[q0][q][pos];
			NewProbs[q] = p;
			NewStates[i] = q;
		}
	}
	else{
		nstates = 0;
		for(k = 0; k < npredstates; k++){
			q1 = PredStates[k];
			s = MainData::ND_HHMTrans[q1][pos].size();
			for(i = 0; i < s; i++){
				q = MainData::ND_HHMTrans[q1][pos][i];
				p = MainData::ND_HHMProbs[q1][q][pos];
				if(NewProbs[q] == 0){
					NewStates[nstates] = q;
					nstates ++;
				}
				NewProbs[q] += PredProbs[q1]*p;
			}
			PredProbs[q1] = 0;
			PredStates[k] = 0;
		}
	}

	if(wordsize - 1 != from){
		npredstates = nstates;
		for(k = 0; k < nstates; k++){
			q = NewStates[k]; 
			NewStates[k] = 0;
			PredStates[k] = q;
			PredProbs[q] = NewProbs[q];
			NewProbs[q] = 0;	
		}
		AllTermProbs(q0, PredProbs, PredStates, npredstates, NewProbs, NewStates, word, wordsize, from + 1);
	}
	return;
};

double ND_HHM_Prob::TermProb(std::string word, int q){
	size_t s = word.size();
    double p = TermCondProb(0,word,s,0,q);
	return p;
};


double ND_HHM_Prob::TransitionProb(int q1, int q2){
	double p = 0;
	int i;
	for(i = 0; i < MainData::AlpSize; i++){
		p += MainData::ND_HHMProbs[q1][q2][i];
	}
	return p;
};

vector<int> ND_HHM_Prob::ConsistStates(int q){
	vector<int> vec;
	int i,j;
	for(i = 0; i < H_M_Node::NumAllStates; i++){
		for(j = 0; j < MainData::AlpSize; j++){
			if(MainData::ND_HHMProbs[i][q][j] != 0){
				vec.push_back(i);
				break;
			}
		}
	}
	return vec;
};


/////////////////////////////////////////////////




void Forward(int nseq, int pos, char** seq, double*** F, double*** HMMProbs){
	int  q;
	
	if(pos == 0){
		F[nseq][0][0] = 1;
		for(q = 1; q < H_M_Node::NumAllStates; q++){
			if(nseq == 0)
				F[0][q][0] = 0;
			else
				F[nseq][q][0] = 1;
		}	
	}
	else{
		int symb = MainData::AToi(seq[nseq][pos - 1]);
		
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			int q1;
			for ( q1 = 0; q1 < H_M_Node::NumAllStates; q1++){ 
				F[nseq][q][pos] +=  F[nseq][q1][pos-1]* HMMProbs[q1][q][symb];
			}
		}	
	}
	if(pos < (int)strlen(seq[nseq])) 
		Forward(nseq, pos + 1, seq, F, HMMProbs);
}




void Backward(int nseq, int pos,  char** seq, double*** B, double*** HMMProbs){
	int  q;
	int seqlen = strlen(seq[nseq]);
	if(pos == seqlen){
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			B[nseq][q][seqlen] = 1;
		}	
	}
	else{
		int symb = MainData::AToi(seq[nseq][pos]);
		
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			int q1;
			B[nseq][q][pos] = 0;
			for ( q1 = 0; q1 < H_M_Node::NumAllStates; q1++){ 
				B[nseq][q][pos] +=  HMMProbs[q][q1][symb]*B[nseq][q1][pos+1];
			}
		}	
	}
	if(pos >1) 
		Backward(nseq, pos - 1,  seq, B, HMMProbs);
}

double Baum_Welch(int Nseqs, char** seq, double*** F, double*** B, double*** HMMProbs, double*** NewHMMProbs){
	
	int q, q1, i, n;
	double like = 0; 
	double Prob;
		
	//I. Compute Expect((q1,a,q))

	for(n = 0; n< Nseqs; n++){
		int seqlen = strlen(seq[n]);
				//I. Compute Prob(seq[n])
		Prob = 0;
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			Prob += F[n][q][seqlen];
		}	
		like +=  log(Prob);


		for(q1 = 0; q1 < H_M_Node::NumAllStates; q1++){
			for(q = 0; q < H_M_Node::NumAllStates; q++){
				for(i = 0; i < seqlen; i++){
					int symb = MainData::AToi(seq[n][i]);
					NewHMMProbs[q1][q][symb] += (F[n][q1][i]*HMMProbs[q1][q][symb]*B[n][q][i+1])/Prob;
				}
			}	
		} 
	}

	//II. Compute pi(q1,a,q)
	for(q1 = 0; q1 < H_M_Node::NumAllStates; q1++){
		double Sum = 0;
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				Sum += NewHMMProbs[q1][q][i];
			}
		}

		for(q = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				NewHMMProbs[q1][q][i] = NewHMMProbs[q1][q][i]/Sum;
			}
		}
	} 
	return like;
}


double ND_HHM_Prob::HMM_Estimate_Param(int Nseqs, char** seq, double*** HMMProbs, double*** NewHMMProbs){

	int seqlen;
	int q, q1, i, step, n;
	
	double like, predlike;

	double*** F = new double**[Nseqs]; 
	double*** B = new double**[Nseqs];
	
	for(n = 0; n < Nseqs; n++){
		F[n] = new double*[H_M_Node::NumAllStates]; 
		B[n] = new double*[H_M_Node::NumAllStates];

		seqlen = strlen(seq[n]);
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			F[n][q] = new double[seqlen + 1]; 
			B[n][q] = new double[seqlen + 1]; 
		}
	}

	predlike = -10000001;
	like = -10000000;
	step = 0;
	while((like - predlike > 0.01) && (step < 200)){
		for(n = 0; n < Nseqs; n++){
			seqlen = strlen(seq[n]);
			for(q = 0; q < H_M_Node::NumAllStates; q++){
				for(i = 0; i <= seqlen; i++){
					F[n][q][i] = 0;
					B[n][q][i] = 0;
				}
			}

			Forward(n, 0, seq, F, HMMProbs);
			Backward(n, seqlen, seq, B, HMMProbs);
		}

		predlike = like;
		like = Baum_Welch(Nseqs, seq, F, B, HMMProbs, NewHMMProbs);
		cout<<step<<'\t'<<like<<'\n';

		for(q1 = 0; q1< H_M_Node::NumAllStates; q1++){
			for(q = 0; q < H_M_Node::NumAllStates; q++){
				for(i = 0; i <= MainData::AlpSize; i++){
					HMMProbs[q1][q][i] =  NewHMMProbs[q1][q][i];
					NewHMMProbs[q1][q][i] = 0;
				}
			}
		}
		step ++;
	}

	/*
	for(n = 0; n < Nseqs; n++){
		seqlen = strlen(seq[n]);
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			delete[] F[n][q]; 
			delete[] B[n][q]; 
			F[n][q] = nullptr; 
			B[n][q] = nullptr; 
		}
		delete[] F[n]; 
		delete[] B[n];
		F[n] = nullptr; 
		B[n] = nullptr;
	}
	delete[] F;
	delete[] B;
	F = nullptr;
	B = nullptr;
	*/
	return like;
}

////////////////////////////////////////////////////
double Rand_in_range(double x){
	double r = (double) rand ()/RAND_MAX;
	r = r* x;
	return r;
}


void RandHMM(double*** HMMProbs){
	time_t time1 = time(NULL);
	int time = int(time1);
	srand(time);
	rand(); 
	
	int q,q1, i;
	for(q1  = 0; q1 < H_M_Node::NumAllStates; q1++){
		double Sum = 0;
		for(q  = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				HMMProbs[q1][q][i] = (double) rand ()/RAND_MAX;;
				Sum += HMMProbs[q1][q][i];
			}
		}

		for(q  = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				HMMProbs[q1][q][i] = HMMProbs[q1][q][i]/Sum;
			}
		}
	}

}


double  ND_HHM_Prob::Find_Best_Model(int NSeqs, char** seq, double*** &BestHMMProbs){
	double*** HMMProbs;
	double*** NewHMMProbs;
	double BestLike = -10000000;

	HMMProbs = new double**[H_M_Node::NumAllStates];
	NewHMMProbs = new double**[H_M_Node::NumAllStates];
	BestHMMProbs = new double**[H_M_Node::NumAllStates];

	int q,q1,ind;
	for(q1 = 0; q1 < H_M_Node::NumAllStates; q1++){
		NewHMMProbs[q1] = new double*[H_M_Node::NumAllStates];
		HMMProbs[q1] = new double*[H_M_Node::NumAllStates];
		BestHMMProbs[q1] = new double*[H_M_Node::NumAllStates];
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			NewHMMProbs[q1][q] = new double[MainData::AlpSize];
			HMMProbs[q1][q] = new double[MainData::AlpSize];
			BestHMMProbs[q1][q] = new double[MainData::AlpSize];
			for(ind = 0; ind < MainData::AlpSize; ind ++){
				NewHMMProbs[q1][q][ind] = 0; 
				HMMProbs[q1][q][ind] = 0; 
				BestHMMProbs[q1][q][ind] = 0; 
			}
		}
	}
	
	int NModels = 5;
	int i;
	for(i  = 0; i < NModels; i++){
		RandHMM(HMMProbs);
		cout<<"Model: "<<i<<'\n';
		double like = ND_HHM_Prob::HMM_Estimate_Param(NSeqs, seq,  HMMProbs, NewHMMProbs);
		if(like > BestLike){
			BestLike = like;
			for(q1 = 0; q1 < H_M_Node::NumAllStates; q1++){
				for(q = 0; q < H_M_Node::NumAllStates; q++){
					for(ind = 0; ind < MainData::AlpSize; ind ++){
						BestHMMProbs[q1][q][ind] = HMMProbs[q1][q][ind]; 
					}
				}
			}
		}
	}



	for(q1  = 0; q1 < H_M_Node::NumAllStates; q1++){
		double Sum = 0;
		for(q  = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				if(BestHMMProbs[q1][q][i] < 0.00000001)
					BestHMMProbs[q1][q][i] = 0;
				Sum += BestHMMProbs[q1][q][i];
			}
		}

		for(q  = 0; q < H_M_Node::NumAllStates; q++){
			for(i = 0; i < MainData::AlpSize; i++){
				BestHMMProbs[q1][q][i] = BestHMMProbs[q1][q][i]/Sum;
			}
		}
	}

	
//	cout<<-2<<'\n';
	//cout<<"A C G T \n";
	//cout<<H_M_Node::NumAllStates<<'\n';
	cout<<"HMM\n";
	for(q1 = 0; q1< H_M_Node::NumAllStates; q1++){
		cout<<"-."<<q1<<'\n';
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			for(ind = 0; ind < MainData::AlpSize; ind++){
				cout<<BestHMMProbs[q1][q][ind]<<'\t';
			}
			cout<<"\n";
		}
	}

	/*
	for(q1 = 0; q1 < H_M_Node::NumAllStates; q1++){
		for(q = 0; q < H_M_Node::NumAllStates; q++){
			delete[] NewHMMProbs[q1][q];
			delete[] HMMProbs[q1][q];
			NewHMMProbs[q1][q] = nullptr;
			HMMProbs[q1][q] = nullptr;
		}
		delete[] NewHMMProbs[q1];
		delete[] HMMProbs[q1];
		NewHMMProbs[q1] = nullptr;
		HMMProbs[q1] = nullptr;
	}		
	delete[] HMMProbs;
	delete[] NewHMMProbs;
	HMMProbs = nullptr;
	NewHMMProbs = nullptr;
	*/
	return BestLike;
};

void ND_HHM_Prob::GenRandHMM1(double d, double*** HMMProbs){
	
	vector<double> vec;
	/*vec.push_back(0.1);
	vec.push_back(0.2); 
	vec.push_back(0.3);
	vec.push_back(0.4);
	vec.push_back(0.5);
	vec.push_back(0.7);
	vec.push_back(1);
	int v_s = 7;
	*/
	vec.push_back(1);
	vec.push_back(0.25);
	//vec.push_back(0.5); 
	//vec.push_back(0.75);
	//vec.push_back(1);
	int v_s = 2;

	int nstates,ind2;
	int q,q1, i;
	int maxnstates = 26;
	int nserias = 51;

	double dtrans = 1;

	double y = dtrans*MainData::AlpSize;
	int ntrans = (int) y; 

	HMMProbs = new double**[maxnstates];

	int ind;
	for(q1 = 0; q1 < maxnstates; q1++){
		HMMProbs[q1] = new double*[maxnstates];
		for(q = 0; q < maxnstates; q++){
			HMMProbs[q1][q] = new double[MainData::AlpSize];
			for(ind = 0; ind < MainData::AlpSize; ind ++){
				HMMProbs[q1][q][ind] = 0; 
			}
		}
	}

	int* pairs = Malloc<int>(maxnstates);
	int* trans = Malloc<int>(MainData::AlpSize);
	int ns = 0;

	time_t time1 = time(NULL);
	int time = int(time1);
	srand(time);
	rand();


	for(ns = 0; ns < nserias; ns++){
		for(nstates = 2; nstates<maxnstates; nstates++){
			for(ind2 = 0; ind2 < v_s; ind2++){
				 

				d = vec[ind2];
				for(q1 = 0; q1< nstates; q1++){
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i]= 0;
						}
					}
				}


				double x = d*nstates;
				int npairs = (int) x;

				if(npairs < x){
					npairs ++;
				}
				
				for(q1  = 0; q1 < nstates; q1++){
			
					memset(pairs, 0x00, maxnstates* sizeof(int));

					int n = 0;
					while(n < npairs){
						x = (double) rand () /RAND_MAX;
						x = x*nstates;
						int pos = (int) x;
						if(pos != nstates){
							if(pairs[pos] == 0){
								pairs[pos] = 1;
								n++;
							}
						}
					}
		
					double Sum = 0;
					for(q  = 0; q < nstates; q++){
						
						if(pairs[q] != 0){
						
							memset(trans, 0x00, MainData::AlpSize* sizeof(int));
							n = 0;
							while(n < ntrans){
								x = (double) rand () /RAND_MAX;
								x = x*MainData::AlpSize;
								int pos = (int) x;
								if(pos != MainData::AlpSize){
									if(trans[pos] == 0){
										trans[pos] = 1;
										n++;
									}
								}
							}

							for(i = 0; i < MainData::AlpSize; i++){
								HMMProbs[q1][q][i] = 0;
								if(trans[i] != 0){
									HMMProbs[q1][q][i] = (double) rand ()/RAND_MAX;
									Sum += HMMProbs[q1][q][i];
								}
							}
						}
					}
					for(q  = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i] = HMMProbs[q1][q][i]/Sum;
						}
					}
				}


				//open of out file
				std::stringstream ststream;
				ststream << nstates;
				std::string st = ststream.str();
				std::stringstream st1stream;
				st1stream << ns << "/";
				std::string st1 = st1stream.str();
				if(d == 0.1){
					st += "_10";
				}
				if(d == 0.2){
					st += "_20";
				}
				if(d == 0.3){
					st += "_30";
				}
					if(d == 0.4){
					st += "_40";
				}
				if(d == 0.5){
					st += "_50";
				}
				if(d == 0.7){
					st += "_70";
				}
				if(d == 1){
					st += "_100";
				}
				if(d == 0.25){
					st += "_25";
				}
				if(d == 0.75){
					st += "_75";
				}
				

				std::string fname = "Alps/Alp" + st1 + st + ".txt";
				ofstream outalp;
				outalp.open(fname.c_str());
			
				outalp<<"-2"<<'\n';
				outalp<<"A C G T\n";
				outalp<<nstates<<"\n\n";
				for(q1 = 0; q1< nstates; q1++){
					outalp<<'-'<<q1<<'.'<<'\n';
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){
							if(HMMProbs[q1][q][i] < 0.000000001){
								HMMProbs[q1][q][i] = 0;
							}
							outalp.setf(ios_base::fixed);
							outalp<<std::setprecision(10)<<HMMProbs[q1][q][i]<<'\t';
							//outalp<<HMMProbs[q1][q][i]<<'\t';
						}
						outalp<<"\n";
					}
				}
				outalp.close();
			}
		}	
	}
}

void ND_HHM_Prob::GenRandHMM2(double d, double*** HMMProbs){
	
	vector<double> vec;
	

	int nstates;
	int q,q1, i;
	int maxnstates = 26;
	int nserias = 51;
	int npairs = 1;
	double dtrans = 1;
	double x;

	double y = dtrans*MainData::AlpSize;
	int ntrans = (int) y; 

	HMMProbs = new double**[maxnstates];

	int ind;
	for(q1 = 0; q1 < maxnstates; q1++){
		HMMProbs[q1] = new double*[maxnstates];
		for(q = 0; q < maxnstates; q++){
			HMMProbs[q1][q] = new double[MainData::AlpSize];
			for(ind = 0; ind < MainData::AlpSize; ind ++){
				HMMProbs[q1][q][ind] = 0; 
			}
		}
	}

	int* pairs = Malloc<int>(maxnstates);
	int* trans = Malloc<int>(MainData::AlpSize);
	int ns = 0;
	time_t time1 = time(NULL);
	int time = int(time1);
	srand(time);
	rand(); 

	for(ns = 0; ns < nserias; ns++){
		for(nstates = 2; nstates<maxnstates; nstates++){
			
				for(q1 = 0; q1< nstates; q1++){
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i]= 0;
						}
					}
				}


				for(q1  = 0; q1 < nstates; q1++){
			
					memset(pairs, 0x00, maxnstates* sizeof(int));
					
					int n = 0;
					while(n < npairs){
						x = (double) rand () /RAND_MAX;
						x = x*nstates;
						int pos = (int) x;
						if(pos != nstates){
							if(pairs[pos] == 0){
								pairs[pos] = 1;
								n++;
							}
						}
					}
		
					double Sum = 0;
					for(q  = 0; q < nstates; q++){
						
						if(pairs[q] != 0){
						
							memset(trans, 0x00, MainData::AlpSize* sizeof(int));
							n = 0;
							while(n < ntrans){
								x = (double) rand () /RAND_MAX;
								x = x*MainData::AlpSize;
								int pos = (int) x;
								if(pos != MainData::AlpSize){
									if(trans[pos] == 0){
										trans[pos] = 1;
										n++;
									}
								}
							}

							for(i = 0; i < MainData::AlpSize; i++){
								HMMProbs[q1][q][i] = 0;
								if(trans[i] != 0){
									HMMProbs[q1][q][i] = (double) rand ()/RAND_MAX;
									Sum += HMMProbs[q1][q][i];
								}
							}
						}
					}
					for(q  = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i] = HMMProbs[q1][q][i]/Sum;
						}
					}
				}


				//open of out file
				std::stringstream ststream;
				ststream << nstates;
				std::stringstream st1stream;
				st1stream << ns << "/";
				std::stringstream st2stream;
				st2stream << "_" << npairs;
				std::string st = ststream.str();
				std::string st1 = st1stream.str();
				std::string st2 = st2stream.str();
				if(d == 0.1){
					st += "_10";
				}
				if(d == 0.2){
					st += "_20";
				}
				if(d == 0.3){
					st += "_30";
				}
					if(d == 0.4){
					st += "_40";
				}
				if(d == 0.5){
					st += "_50";
				}
				if(d == 0.7){
					st += "_70";
				}
				if(d == 1){
					st += "_100";
				}
				if(d == 0.25){
					st += "_25";
				}
				if(d == 0.75){
					st += "_75";
				}
				

				std::string fname = "Alps/Alp" + st1 + st + st2 + ".txt";
				ofstream outalp;
				outalp.open(fname.c_str());
			
				outalp<<"-2"<<'\n';
				outalp<<"A C G T\n";
				outalp<<nstates<<"\n\n";
				for(q1 = 0; q1< nstates; q1++){
					outalp<<'-'<<q1<<'.'<<'\n';
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){
							if(HMMProbs[q1][q][i] < 0.000000001){
								HMMProbs[q1][q][i] = 0;
							}
							outalp.setf(ios_base::fixed);
							outalp<<std::setprecision(10)<<HMMProbs[q1][q][i]<<'\t';
							//outalp<<HMMProbs[q1][q][i]<<'\t';
						}
						outalp<<"\n";
					}
				}
				outalp.close();
		}	
	}
}
/*
void ND_HHM_Prob::GenRandHMM2(double d, double*** HMMProbs){
	
	vector<double> vec;
	//vec.push_back(0.1);
	vec.push_back(0.25);
	vec.push_back(0.5); 
	vec.push_back(0.75);
	vec.push_back(1);
	int v_s = 4;

	int nstates,ind2;
	int q,q1, i;
	int maxnstates = 26;
	int nserias = 21;

	double dtrans = 0.5;

	double y = dtrans*MainData::AlpSize;
	int ntrans = (int) y; 

	HMMProbs = new double**[maxnstates];

	int ind;
	for(q1 = 0; q1 < maxnstates; q1++){
		HMMProbs[q1] = new double*[maxnstates];
		for(q = 0; q < maxnstates; q++){
			HMMProbs[q1][q] = new double[MainData::AlpSize];
			for(ind = 0; ind < MainData::AlpSize; ind ++){
				HMMProbs[q1][q][ind] = 0; 
			}
		}
	}

	int* pairs = Malloc<int>(maxnstates);
	int* trans = Malloc<int>(MainData::AlpSize);
	int ns = 0;
	for(ns = 0; ns < nserias; ns++){
		for(nstates = 2; nstates<maxnstates; nstates++){
				time_t time1 = time(NULL);
				int time = int(time1);
				srand(time);
				rand(); 

				for(q1 = 0; q1< nstates; q1++){
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i]= 0;
						}
					}
				}


				double x;
				double Sum;
				for(q1  = 0; q1 < nstates; q1++){
					Sum = 0;
					for(i = 0; i < MainData::AlpSize; i++){
						x = (double) rand () /RAND_MAX;
						x = x*nstates;
						int state = (int) x;
						if(state == nstates){
							state --;
						}


						HMMProbs[q1][state][i] = 0;
						HMMProbs[q1][state][i] = (double) rand ()/RAND_MAX;
						Sum += HMMProbs[q1][state][i];
					}
					for(q  = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i] = HMMProbs[q1][q][i]/Sum;
						}
					}
				}


				//open of out file
				char buf[3];
				std::string st;
				st += itoa(nstates,buf,10);
				char buf1[3];
				std::string st1;
				st1 += itoa(ns,buf1,10);
				st1 += "/";
				

				std::string fname = "Alps/Alp" + st1 + st + ".txt";
				ofstream outalp;
				outalp.open(fname.c_str());
			
				outalp<<"-2"<<'\n';
				outalp<<"A C G T\n";
				outalp<<nstates<<"\n\n";
				for(q1 = 0; q1< nstates; q1++){
					outalp<<'-'<<q1<<'.'<<'\n';
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){
							if(HMMProbs[q1][q][i] < 0.000000001){
								HMMProbs[q1][q][i] = 0;
							}
							outalp.setf(ios_base::fixed);
							outalp<<std::setprecision(10)<<HMMProbs[q1][q][i]<<'\t';
							//outalp<<HMMProbs[q1][q][i]<<'\t';
						}
						outalp<<"\n";
					}
				}
				outalp.close();
		}	
	}
}

*/



void ND_HHM_Prob::GenRandHMM(double d, double*** HMMProbs){
	
	vector<double> vec;
	/*vec.push_back(0.1);
	vec.push_back(0.2); 
	vec.push_back(0.3);
	vec.push_back(0.4);
	vec.push_back(0.5);
	vec.push_back(0.7);
	vec.push_back(1);
	int v_s = 7;
	*/
	vec.push_back(0.1);
	vec.push_back(0.25);
	vec.push_back(0.5); 
	vec.push_back(0.75);
	vec.push_back(1);
	int v_s = 5;

	int nstates,ind2;
	int q,q1, i;
	int maxnstates = 26;
	int nserias = 21;

	HMMProbs = new double**[maxnstates];

	int ind;
	for(q1 = 0; q1 < maxnstates; q1++){
		HMMProbs[q1] = new double*[maxnstates];
		for(q = 0; q < maxnstates; q++){
			HMMProbs[q1][q] = new double[MainData::AlpSize];
			for(ind = 0; ind < MainData::AlpSize; ind ++){
				HMMProbs[q1][q][ind] = 0; 
			}
		}
	}

	int* positions = Malloc<int>(maxnstates*MainData::AlpSize);
	int ns = 0;
	for(ns = 0; ns < nserias; ns++){
		for(nstates = 2; nstates<maxnstates; nstates++){
			for(ind2 = 0; ind2 < v_s; ind2++){
				time_t time1 = time(NULL);
				int time = int(time1);
				srand(time);
				rand(); 

				d = vec[ind2];
				for(q1 = 0; q1< nstates; q1++){
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i]= 0;
						}
					}
				}

				double x = d*nstates*MainData::AlpSize;
				int npos = (int) x;

	
				for(q1  = 0; q1 < nstates; q1++){
			
					memset(positions, 0x00, maxnstates*MainData::AlpSize* sizeof(int));

					int n = 0;
					while(n < npos){
						x = (double) rand () /RAND_MAX;
						x = x*nstates*MainData::AlpSize;
						int pos = (int) x;
						if(positions[pos] == 0){
							positions[pos] = 1;
							n++;
						}
					}
		
					double Sum = 0;
					for(q  = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){
							HMMProbs[q1][q][i] = 0;
							if(positions[q*MainData::AlpSize + i] != 0){
								HMMProbs[q1][q][i] = (double) rand ()/RAND_MAX;
								Sum += HMMProbs[q1][q][i];
							}
						}
					}

					for(q  = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){	
							HMMProbs[q1][q][i] = HMMProbs[q1][q][i]/Sum;
						}
					}
				}


				//open of out file
				std::stringstream ststream;
				ststream << nstates;
				std::string st = ststream.str();
				std::stringstream st1stream;
				st1stream << ns << "/";
				std::string st1 = st1stream.str();
				if(d == 0.1){
					st += "_10";
				}
				if(d == 0.2){
					st += "_20";
				}
				if(d == 0.3){
					st += "_30";
				}
					if(d == 0.4){
					st += "_40";
				}
				if(d == 0.5){
					st += "_50";
				}
				if(d == 0.7){
					st += "_70";
				}
				if(d == 1){
					st += "_100";
				}
				if(d == 0.25){
					st += "_25";
				}
				if(d == 0.75){
					st += "_75";
				}
				

				std::string fname = "Alps/Alp" + st1 + st + ".txt";
				ofstream outalp;
				outalp.open(fname.c_str());
			
				outalp<<"-2"<<'\n';
				outalp<<"A C G T\n";
				outalp<<nstates<<"\n\n";
				for(q1 = 0; q1< nstates; q1++){
					outalp<<'-'<<q1<<'.'<<'\n';
					for(q = 0; q < nstates; q++){
						for(i = 0; i < MainData::AlpSize; i++){
							if(HMMProbs[q1][q][i] < 0.000000001){
								HMMProbs[q1][q][i] = 0;
							}
							outalp.setf(ios_base::fixed);
							outalp<<std::setprecision(10)<<HMMProbs[q1][q][i]<<'\t';
							//outalp<<HMMProbs[q1][q][i]<<'\t';
						}
						outalp<<"\n";
					}
				}
				outalp.close();
			}
		}	
	}
}


