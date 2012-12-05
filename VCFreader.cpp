/*
 * VCFreader
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "VCFreader.h"


VCFreader::VCFreader(string file,string indexForFile,string chrName,int start,int end,int indelsAhead){
    readAhead=indelsAhead;
    rt =new ReadTabix (file,indexForFile,chrName,start,end);

    needToPopulateQueue=true;
    fullQueue          =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    svcfToReturn=0;

    
    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    tabixMode = true;
    textMode  = false;
}




VCFreader::VCFreader(string file,int indelsAhead=5){
    readAhead=indelsAhead;
    numberOfTimesHasDataWasCalled=0;
    svcfToReturn=0;

    vcfFile.open(file.c_str(), ios::in);    // open the streams
    if (vcfFile.is_open()) {
	//fine
    }else{
	cerr<<"Unable to open the file "<<file<<endl;
	exit(1);

    }

    needToPopulateQueue =true;
    fullQueue           =false;
    endQueue            =false;


    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    tabixMode = false;
    textMode  = true;
}




VCFreader::~VCFreader(){

    if(tabixMode){
	delete rt; //calling the destructor
    }
    if(textMode)
	vcfFile.close();
    queueOfVCFs.clear();

    delete svcfToReturn;
}

void VCFreader::repositionIterator(string chrName,int start,int end){

    if(!tabixMode){
	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
	exit(1);	
    }


    //re-initializing variables
    needToPopulateQueue=true;
    fullQueue          =false;
    endQueue            =false;
    numberOfTimesHasDataWasCalled=0;

    delete svcfToReturn;
    svcfToReturn=0;

    
    indexInQueueOfIndels=-1;
    indexOfLastIndel=0;
    previouslyFoundIndel=false;

    queueOfVCFs.clear();

    rt->repositionIterator(chrName,start,end);

}

bool VCFreader::getNextLine(){
    if(tabixMode){
	return rt->readLine(currentline);
    }

    if(textMode){
	while(1){
	    bool flag=getline(vcfFile,currentline);
	    if(!flag)
		return false;	   
	    if(currentline.length() > 0 && currentline[0] != '#')
		return true;
	}
    }
}


bool VCFreader::hasData(){

    numberOfTimesHasDataWasCalled++;

    
    //if first call and queue empty, populate
    if(needToPopulateQueue){
	bool loop=true;
	int indexQueue=0;
	while(loop){
	    if(getNextLine()){
		SimpleVCF svcf (currentline);
		if(queueOfVCFs.size() != 0 ){
		    flagCpG( &(queueOfVCFs.back()),&svcf);
		}
		//cout<<"Adding "<<svcf<<endl;
		queueOfVCFs.push_back(svcf);

		if(svcf.containsIndel()){
		    if(indexInQueueOfIndels == -1)
			indexInQueueOfIndels=indexQueue;		   
		}
		indexQueue++;
		if(queueOfVCFs.size() == (readAhead+1)){
		    loop=false;
		}
	    }else{
		loop=false;
	    }
	}
	

	if(queueOfVCFs.size() == (readAhead+1)){
	    fullQueue=true;
	}else{
	    endQueue=true;
	}

	needToPopulateQueue=false;
    }

    //if subsequent call, and queue full
    if(fullQueue){
	if(getNextLine()){
	    SimpleVCF svcf (currentline);
	    if(queueOfVCFs.size() != 0 ){
		flagCpG( &(queueOfVCFs.back()),&svcf);
	    }	    
	    queueOfVCFs.push_back(svcf);
	    
	    if(svcf.containsIndel()){	
		if(indexInQueueOfIndels == -1)
		  indexInQueueOfIndels=queueOfVCFs.size()-1;
	    }
	    
	}else{
	    fullQueue=false;
	    endQueue=true;
	}

    }

    //if final calls and queue not max size
    if(endQueue){
	//nothing to do
    }


    return !(queueOfVCFs.empty());
}


inline bool VCFreader::flagCpG(SimpleVCF * previous,SimpleVCF * current){ //pass by address
    if( ( (previous->getPosition()+1) == current->getPosition())  &&   //one position behind
	(  previous->getChr()         == current->getChr()    )   &&   //on same chr
	(previous->hasAtLeastOneC()   && current->hasAtLeastOneG()) ){  //previous has at least one C, current has at least one G
	previous->setCpg(true);
	current->setCpg(true);
    }
}



SimpleVCF * VCFreader::getData(){
// auto_ptr<AlleleInfo> VCFreader::getData(){

    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;

    //delete the previous data
    delete svcfToReturn;
	//}

    //SimpleVCF svcf =queueOfVCFs.front();

    svcfToReturn =new SimpleVCF(queueOfVCFs.front());
    // cout<<&(queueOfVCFs.front())<<endl;
    // cout<<svcfToReturn<<endl;
    //delete svcfToReturn;
    queueOfVCFs.pop_front();



    //look ahead   
    if(indexInQueueOfIndels == 0){//last element that is to be return is an indel
	svcfToReturn->setCloseIndel(true);
	indexInQueueOfIndels = -1;
    }else{
	if(indexInQueueOfIndels != 0){//elements in the list that are indels, need to check them
	    int indexInList=0;	    
	    list<SimpleVCF>::iterator it;
	    for ( it=queueOfVCFs.begin() ; it != queueOfVCFs.end(); it++ ){
		if(indexInList<=indexInQueueOfIndels){
		    if( it->containsIndel()){

			if( ((it->getPosition() - svcfToReturn->getPosition() ) <= readAhead) &&
			    (it->getChr() == svcfToReturn->getChr()) ){
			    svcfToReturn->setCloseIndel(true);			    
			}

		    }
		}else{
		    break;
		}
		
		indexInList++;
	    }

	}
	indexInQueueOfIndels=max(indexInQueueOfIndels-1,-1);
	
    }


    //look behind for indels
    if(previouslyFoundIndel){
	 if( ((svcfToReturn->getPosition() - indexOfLastIndel) <= readAhead) &&
	     (chrOfLastIndel == svcfToReturn->getChr()) ){
	     svcfToReturn->setCloseIndel(true);
	 }
     }

     

    if(svcfToReturn->containsIndel()){
	previouslyFoundIndel=true;
	indexOfLastIndel    =svcfToReturn->getPosition();
	chrOfLastIndel      =svcfToReturn->getChr();
    }

    return svcfToReturn;
}
