/*
 * SetVCFFilters
 * Date: Jan-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SetVCFFilters.h"

SetVCFFilters::SetVCFFilters(int    minGQcutoff         ,
			     int    minMQcutoff         ,
			     double minMapabilitycutoff ,
			     bool   filterIndelProx     ,
			     bool   repeatMasking       ,
			     bool   systemError         ,    
			     int    minCovcutoff        ,
			     int    maxCovcutoff        ){
    

    this->minGQcutoff          = minGQcutoff         ;
    this->minMQcutoff          = minMQcutoff     ;
    this->minMapabilitycutoff  = minMapabilitycutoff;
    //    this->bpForIndels          = bpForIndels;
    this->filterIndelProx      = filterIndelProx;
    this->repeatMasking        = repeatMasking;
    this->systemError          = systemError;

    if(minCovcutoff < 0 ){
	cerr<<"SetVCFFilters: Cannot have a negative minimum cutoff"<<endl;
	exit(1);
    }

    if(maxCovcutoff < 0 ){
	cerr<<"SetVCFFilters: Cannot have a negative maximum cutoff"<<endl;
	exit(1);
    }

    this->minCovcutoff         = minCovcutoff;
    this->maxCovcutoff         = maxCovcutoff;
    this->name ="";
}

SetVCFFilters::~SetVCFFilters(){

}

void SetVCFFilters::setName(string name) {
    this->name = name;
}

string SetVCFFilters::getName() const{
    return name;
}

int SetVCFFilters::getMinGQcutoff() const {
    return minGQcutoff;
}

int SetVCFFilters::getMinMQcutoff() const {
    return minMQcutoff;
}

double SetVCFFilters::getMinMapabilitycutoff() const {
    return minMapabilitycutoff;
}

// int SetVCFFilters::getBpForIndels() const {
//     return bpForIndels;
// }

bool SetVCFFilters::getFilterIndelProx() const {
    return filterIndelProx;
}

bool SetVCFFilters::getRepeatMasking() const {
    return repeatMasking;
}

bool SetVCFFilters::getSystemError() const {
    return systemError;
}

int  SetVCFFilters::getMinCovcutoff() const {
    return minCovcutoff;
}

int  SetVCFFilters::getMaxCovcutoff() const {
    return maxCovcutoff;
}
