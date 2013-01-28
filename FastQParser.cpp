/*
 * FastQParser
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FastQParser.h"

FastQParser::FastQParser(string file){
    //cout<<file<<endl;
    fastqFile.open(file.c_str(), ios::in);    // open the streams

    if (fastqFile.is_open()) {
	//fine
    }else{
	cerr<<"Unable to open the file "<<file<<endl;
	exit(1);

    }
    numberOfTimesHasDataWasCalled=0;
    initPointer=false;
}

FastQParser::~FastQParser(){
    fastqFile.close();
    // delete toreturn;
}


bool       FastQParser::hasData(){

    if(!initPointer){
	initPointer=true;
    }else{
	// cout<<"hasData"<<endl;
	delete(toreturn);
    }

    numberOfTimesHasDataWasCalled++;
    bool flag;

    flag=getline(fastqFile,currentline1);

    if(!flag){
	return false;
    }else{	
	if(currentline1.length() > 0 && currentline1[0] != '@'){
	    cerr<<"Parse error with line "<<currentline1<<endl;
	    exit(1);
	}
    }
    // cout<<"1 "<<currentline1<<endl;
    flag=getline(fastqFile,currentline2);
    // cout<<"2 "<<currentline2<<endl;
    if(!flag){
	cerr<<"Parse error with line "<<currentline2<<endl;
	exit(1);
    }

    flag=getline(fastqFile,currentline3);
    // cout<<"3 "<<currentline3<<endl;
    if(!flag){
	cerr<<"Parse error with line "<<currentline3<<endl;
	exit(1);
    }else{	
	if(currentline3.length() > 0 && currentline3[0] != '+'){
	    cerr<<"Parse error with line "<<currentline3<<endl;
	    exit(1);
	}
    }

    flag=getline(fastqFile,currentline4);
    // cout<<"4 "<<currentline4<<endl;
    if(!flag){
	cerr<<"Parse error with line "<<currentline4<<endl;
	exit(1);
    }

    // cout<<"new "<<endl;
    toreturn=new FastQObj(&currentline1,&currentline2,&currentline4);
    return true;
}

FastQObj  *     FastQParser::getData(){

    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;

    return toreturn;
}
