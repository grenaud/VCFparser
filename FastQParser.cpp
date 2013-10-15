/*
 * FastQParser
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FastQParser.h"

FastQParser::FastQParser(string file,bool isFasta){
    this->isFasta = isFasta;
    
    fastqFile.open(file.c_str(), ios::in);    // open the streams

    if (fastqFile.good()) {
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
	delete(toreturn);
    }

    numberOfTimesHasDataWasCalled++;
    bool flag;


    if(isFasta){
	flag=getline(fastqFile,currentline1);
	if(!flag){
	    return false;
	}else{	
	    if(currentline1.length() > 0 && currentline1[0] != '>'){
		cerr<<"Parse error#1 with line "<<currentline1<<endl;
		exit(1);
	    }
	}

	string * id = new string(currentline1);
	string * seq = new string();

	while(flag){	
	    char peekChar=fastqFile.peek();
	    if(peekChar ==  EOF){
		toreturn=new FastQObj(id,seq,0,isFasta);		  
		return true;
	    }

	    if(peekChar ==  '>'){//beginning of next record
		toreturn=new FastQObj(id,seq,0,isFasta);		      
		return true;
	    }


	    flag=getline(fastqFile,currentline2);

	    if(!flag){
		cerr<<"Parse error#2 with line "<<currentline1<<endl;
		exit(1);
	    }else{	
		seq->append(currentline2);		    
	    }
	}

	

    }else{
	flag=getline(fastqFile,currentline1);

	if(!flag){
	    return false;
	}else{	
	    if(currentline1.length() > 0 && currentline1[0] != '@'){
		cerr<<"Parse error#3 with line "<<currentline1<<endl;
		exit(1);
	    }
	}

	flag=getline(fastqFile,currentline2);

	if(!flag){
	    cerr<<"Parse error#4 with line "<<currentline2<<endl;
	    exit(1);
	}

	flag=getline(fastqFile,currentline3);

	if(!flag){
	    cerr<<"Parse error#5 with line "<<currentline3<<endl;
	    exit(1);
	}else{	
	    if(currentline3.length() > 0 && currentline3[0] != '+'){
		cerr<<"Parse error#6 with line "<<currentline3<<endl;
		exit(1);
	    }
	}

	flag=getline(fastqFile,currentline4);

	if(!flag){
	    cerr<<"Parse error#7 with line "<<currentline4<<endl;
	    exit(1);
	}


	toreturn=new FastQObj(&currentline1,&currentline2,&currentline4,isFasta);
    }
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
