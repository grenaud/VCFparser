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

public:
    FastQObj(string * _id,string * _seq,string * _qual);
    ~FastQObj();
    
    string * getSeq() const;
    string *  getID() const;
    string * getQual() const;
    
    friend ostream& operator<<(ostream& str, FastQObj const& fqo){
	str<<*(fqo.id)<<endl<<*(fqo.seq)<<endl<<"+"<<endl<<*(fqo.qual);
	return str;
    }

};
#endif
