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

#include "FastQObj.h"


using namespace std;

class FastQParser{
private:
    ifstream fastqFile; 
    int numberOfTimesHasDataWasCalled;
    string currentline1;
    string currentline2;
    string currentline3;
    string currentline4;

    FastQObj * toreturn;

public:
    FastQParser(string file);
    ~FastQParser();

    bool       hasData();
    FastQObj * getData();
    bool initPointer;
};
#endif
