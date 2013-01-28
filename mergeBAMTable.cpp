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


    for(int i=1;i<=(argc-1);i++){
	//cout<<argv[i]<<endl;
	// if(!isFile(argv[i])){
	//     cerr << "File "<<argv[1] <<" is not a file"<<endl;
	//     return 1;       
	// }
	BAMTABLEreader * bt =  new BAMTABLEreader(argv[i]);
	btrvec.push_back( bt  );
    }


    BAMTableObj   * btObjvec [btrvec.size()];
    bool    haveData [btrvec.size()];
    
    //find min coordinate
    int minCoord=-1;
    for(unsigned int i=0;i<btrvec.size();i++){    
	haveData[i]=btrvec[i]->hasData();
	if(!haveData[i]){
	    cout<<"File "<<argv[i+1]<<" does not have data"<<endl;
	    return 1;
	}
	btObjvec[i]=btrvec[i]->getData();
	if(i==0){
	    minCoord=btObjvec[i]->getPosition();
	}else{
	    if(int(btObjvec[i]->getPosition())<minCoord ){
		minCoord=int(btObjvec[i]->getPosition());
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
		    }

		}
	    }
	}
	

	if(!hasAtLeastOneWithData)
	    break;
	if(entered)
	    cout<<concat<<endl;
	// return 1;

	currentCoord++;

    }


    
    
    return 0;
}

