#include <iostream>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<armadillo>
#include<cmath>
#include <cstdlib>


using namespace std;

struct solution {
             float x;
             float v;
            };
            
//Count rows
int flcount(ifstream& in){
	                  string s;
                          int k=0;

                          getline(in,s);
                          while (!in.eof()) {
                                             k++;
                                             getline(in,s);
                                            }

                          in.clear();
                          in.seekg(0);
			  return k;
                         }            

int main() {
	    	ifstream in;
        	in.open("solutions_7_10pts.txt");
            
            if ( in.fail() ) {
            cout << "Attention! problems with file opening "<< endl;
            return 1; 
            }
            
            
//Caricamento array e stampa risultati
            int k = flcount(in);
            solution *S = new solution [k];
        
            for(int i=0; i < k; i++){
				     in>>S[i].x;
                      in>>S[i].v;
                                     }
                                     
             cout<<"I dati caricati nell'array sono: "<<endl;
            for(int i=0; i<k; i++){
                                                      cout<<S[i].x<<"\t"<<S[i].v<<endl;
                                                      }
                               cout<<endl;
                                     
return 0;


}


                                     
