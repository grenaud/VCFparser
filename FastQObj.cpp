/*
 * FastQObj
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FastQObj.h"

FastQObj::FastQObj(string * _id,string * _seq,string * _qual,bool _isFasta){

    id   =_id;
    seq  =_seq;
    qual =_qual;
    isFasta =_isFasta;
    if(!isFasta){
	if(seq->length() != qual->length() ){
	    cerr<<"Length of sequence is not the same as the length of the quality"<<endl;
	    exit(1);
	}
    }

}

FastQObj::~FastQObj(){
    // cout<<"dest"<<endl;
    //
    if(isFasta){
	delete(id);
	delete(seq);
    }
    //	delete(qual);
}


void FastQObj::printFastaSeqWithBreaks(ostream & stro) const{
    stro<<*id<<endl;
    for(unsigned int i=0;
	i<seq->size();
	i++){
	if(i%60==0 && i!=0){
	    stro<<endl;
	}

	stro<<seq->at(i);
    }
    stro<<endl;
}

string * FastQObj::getSeq()  const{  return seq;  }
string * FastQObj::getID()   const{  return id;   }

void FastQObj::setID(const string  newID){
    delete(id);
    id=new string(newID);
}

string * FastQObj::getQual() const{  
    if(isFasta){
	cerr<<"Record for id = "<<(*id)<< "is of type fasta" <<endl;
	exit(1);
    }
    return qual;
}
