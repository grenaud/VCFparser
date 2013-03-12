/*
 * mergeBAMTable
 * Date: Oct-26-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "BAMTABLEreader.h"
#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

   // string line;
   // ifstream myFile;
   // string filename = string(argv[1]);
   // myFile.open(filename.c_str(), ios::in);

   // if (myFile.is_open()){
   //   while ( getline (myFile,line)){

   //   }
   //   myFile.close();
   // }else{
   //     cerr << "Unable to open file "<<filename<<endl;
   //     return 1;
   //  }


    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<argv[0] <<" [bam table1] [bam table2] [bam table3] .. "<<"prints to stdout"<<endl;
	return 1;       
    }

    vector<BAMTABLEreader* > btrvec;
    unsigned int total=0;

    for(int i=1;i<=(argc-1);i++){
	//cout<<argv[i]<<endl;
	// if(!isFile(argv[i])){
	//     cerr << "File "<<argv[1] <<" is not a file"<<endl;
	//     return 1;       
	// }
       	    cerr << "Using file "<<argv[1] <<""<<endl;
	BAMTABLEreader * bt =  new BAMTABLEreader(argv[i]);
	btrvec.push_back( bt  );
    }


    BAMTableObj   * btObjvec [btrvec.size()];
    bool    haveData [btrvec.size()];
    string chr;
    //find min coordinate
    int minCoord=-1;
    for(unsigned int i=0;i<btrvec.size();i++){    
	haveData[i]=btrvec[i]->hasData();
	if(!haveData[i]){
	    cerr<<"File "<<argv[i+1]<<" does not have data"<<endl;
	    //return 1;
	    continue;
	}
	btObjvec[i]=btrvec[i]->getData();
	if(i==0){
	    minCoord=btObjvec[i]->getPosition();
	    chr    =btObjvec[i]->getChr();
	}else{
	    if(int(btObjvec[i]->getPosition())<minCoord ){
		minCoord=int(btObjvec[i]->getPosition());
		if(chr   != btObjvec[i]->getChr()){
		    cout<<"File "<<argv[i+1]<<" has a different chr"<<endl;
		    return 1;
		}
	    }
	}
    }
    

    int currentCoord=minCoord;
    while(1){
	BAMTableObj concat;
	bool entered=false;
	bool hasAtLeastOneWithData=false;
	// int todelete=0;

	for(unsigned int i=0;i<btrvec.size();i++){    
	    hasAtLeastOneWithData |= haveData[i];
	    if(haveData[i]){
		if(currentCoord  == int(btObjvec[i]->getPosition()) ){
		    // todelete++;
		    entered=true;
		    concat+=*(btObjvec[i]);
		    haveData[i]=btrvec[i]->hasData();
		    if(haveData[i]){
			btObjvec[i]=btrvec[i]->getData();
			if(chr   != btObjvec[i]->getChr()){
			    cout<<"File "<<argv[i+1]<<" has a different chr"<<endl;
			    return 1;
			}
		    }

		}
	    }
	}
	

	if(!hasAtLeastOneWithData)
	    break;
	if(entered){
	    cout<<concat<<endl;
	    total++;
	}
	// return 1;

	currentCoord++;

    }


    cerr<<"Terminated gracefully, wrote "<<total<<" records "<<endl;
    
    return 0;
}

