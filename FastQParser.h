/*
 * FastQParser
 * Date: Dec-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef FastQParser_h
#define FastQParser_h

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <gzstream.h>

#include "FastQObj.h"


using namespace std;

class FastQParser{
private:
    bool isFasta;
    igzstream fastqFile; 
    int numberOfTimesHasDataWasCalled;
    string currentline1;
    string currentline2;
    string currentline3;
    string currentline4;

    FastQObj * toreturn;

public:
    FastQParser(string file,bool isFasta=false);
    ~FastQParser();

    bool       hasData();
    FastQObj * getData();
    bool initPointer;
};
#endif
