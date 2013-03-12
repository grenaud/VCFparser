/*
 * VCFreader
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "BAMTABLEreader.h"


BAMTABLEreader::BAMTABLEreader(string file,string indexForFile,string chrName,int start,int end){
    readAhead=1;//for CpGs, can be increased
    rt =new ReadTabix (file,indexForFile,chrName,start,end);
    needToPopulateQueue=true;
    fullQueue          =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    btoToReturn=0;

    

    tabixMode = true;
    textMode  = false;
}




BAMTABLEreader::BAMTABLEreader(string file,int indelsAhead){
    readAhead=indelsAhead;
    numberOfTimesHasDataWasCalled=0;
    btoToReturn=0;

    bmtblFile.open(file.c_str(), ios::in);    // open the streams
    if (bmtblFile.good()) {
	//fine
    }else{
	cerr<<"Unable to open the file "<<file<<endl;
	exit(1);
    }

    needToPopulateQueue =true;
    fullQueue           =false;
    endQueue            =false;


    tabixMode = false;
    textMode  = true;
}




BAMTABLEreader::~BAMTABLEreader(){

    if(tabixMode)
	delete rt; //calling the destructor
    if(textMode)
	bmtblFile.close();
    //if( !(btoToReturn) ){
    // cout<<"DESTROYVCF"<<endl;
    queueOfBTOs.clear();
    delete btoToReturn;
    //}
}


void BAMTABLEreader::repositionIterator(string chrName,int start,int end){
    if(!tabixMode){
	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
	exit(1);	
    }

    //re-initializing variables
    needToPopulateQueue=true;
    fullQueue          =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    delete btoToReturn;
    btoToReturn=0;

    queueOfBTOs.clear();

    rt->repositionIterator(chrName,start,end);    
  
}

bool BAMTABLEreader::getNextLine(){
    if(tabixMode){
	return rt->readLine(currentline);
    }

    if(textMode){
	while(1){
	    bool flag=getline(bmtblFile,currentline);

	    if(!flag)
		return false;	   
	    if(currentline.length() > 0 && currentline[0] != '#')
		return true;
	}
    }

    cerr<<"Invalid state in BAMTABLEreader::getNextLine()"<<endl;
    exit(1);
    return false;
}



inline void BAMTABLEreader::flagCpG(BAMTableObj * previous,BAMTableObj * current){ //pass by address
    if( ( (previous->getPosition()+1) == current->getPosition())  &&   //one position behind
	(  previous->getChr()         == current->getChr()    )   &&   //on same chr
	(previous->hasAtLeastOneC()   && current->hasAtLeastOneG()) ){  //previous has at least one C, current has at least one G
	previous->setCpg(true);
	current->setCpg(true);
    }    
}


bool BAMTABLEreader::hasData(){

    numberOfTimesHasDataWasCalled++;

    
    //if first call and queue empty, populate
    if(needToPopulateQueue){
	bool loop=true;

	while(loop){
	    if(getNextLine()){
		BAMTableObj svcf (currentline);
		//detect CpGs
		if(queueOfBTOs.size() != 0 ){
		    flagCpG( &(queueOfBTOs.back()),&svcf);
		}

		queueOfBTOs.push_back(svcf);

		if(queueOfBTOs.size() == (readAhead+1)){
		    loop=false;
		}
	    }else{
		loop=false;
	    }
	}
	

	if(queueOfBTOs.size() == (readAhead+1)){
	    fullQueue=true;
	}else{
	    endQueue=true;  
	}

	needToPopulateQueue=false;
    }

    //if subsequent call, and queue full
    if(fullQueue){
	if(getNextLine()){
	    BAMTableObj svcf (currentline);	    
	    //detect CpGs
	    if(queueOfBTOs.size() != 0 ){
		flagCpG( &(queueOfBTOs.back()),&svcf);
	    }

	    queueOfBTOs.push_back(svcf);
	    
	
	    
	}else{
	    fullQueue=false;
	    endQueue=true;
	}

    }

    //if final calls and queue not max size
    if(endQueue){
	//nothing to do
    }


    return !(queueOfBTOs.empty());
}



BAMTableObj * BAMTABLEreader::getData(){

    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;

    //delete the previous data
    delete btoToReturn;

    btoToReturn =new BAMTableObj(queueOfBTOs.front());
    queueOfBTOs.pop_front();


    return btoToReturn;
}
