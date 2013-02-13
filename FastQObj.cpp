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


string * FastQObj::getSeq()  const{  return seq;  }
string * FastQObj::getID()   const{  return id;   }
string * FastQObj::getQual() const{  
    if(isFasta){
	cerr<<"Record for id = "<<(*id)<< "is of type fasta" <<endl;
	exit(1);
    }
    return qual;
}
