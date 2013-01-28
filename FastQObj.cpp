/*
 * FastQObj
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FastQObj.h"

FastQObj::FastQObj(string * _id,string * _seq,string * _qual){
    if(_seq->length() != _qual->length() ){
	cerr<<"Length of sequence is not the same as the length of the quality"<<endl;
	exit(1);
    }

    id   =_id;
    seq  =_seq;
    qual =_qual;
    
}

FastQObj::~FastQObj(){
    // cout<<"dest"<<endl;
    // delete(id);
    // delete(seq);
    // delete(qual);
}


string * FastQObj::getSeq()  const{  return seq;  }
string * FastQObj::getID()   const{  return id;   }
string * FastQObj::getQual() const{  return qual; }
