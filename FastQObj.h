/*
 * FastQObj
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef FastQObj_h
#define FastQObj_h

#include <string>
#include <iostream>
#include <cstdlib>

using namespace std;

class FastQObj{
private:
    string * id;
    string * seq;
    string * qual;
    bool isFasta;

public:
    FastQObj(string * _id,string * _seq,string * _qual,bool _isFasta=false);
    ~FastQObj();
    
    string * getSeq() const;
    string *  getID() const;
    void setID(const string  newID);

    string * getQual() const;
    void printFastaSeqWithBreaks(ostream & stro) const;

    friend ostream& operator<<(ostream& str, FastQObj const& fqo){
	if(fqo.isFasta)
	    str<<*(fqo.id)<<endl<<*(fqo.seq);
	else
	    str<<*(fqo.id)<<endl<<*(fqo.seq)<<endl<<"+"<<endl<<*(fqo.qual);
	return str;
    }

};
#endif
