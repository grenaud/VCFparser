/*
 * ReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "ReadTabix.h"

// using namespace std;

ReadTabix::ReadTabix(string file,string indexForFile,string chrName,int start,int end){
    if(!isFile(file)){
	cerr<<"File "<<file<<" does not exist"<<endl;
	exit(1);
    }

    if(!isFile(indexForFile)){
	cerr<<"File "<<file<<" does not exist"<<endl;
	exit(1);
    }


    fpTab=ti_open(file.c_str(),indexForFile.c_str());	       

    // -1 is substracted from the start because, for some reason that is unknown to me
    // Heng Li does that in index.c in int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
    // to the get the coordinates right
    // I do this for consistency with the command line program   
    if(start>0)
	start--;

    iteratorTab=ti_query(fpTab,chrName.c_str(),start,end); 
}

ReadTabix::~ReadTabix(){
    ti_iter_destroy(iteratorTab);
    ti_close(fpTab);
}



void ReadTabix::repositionIterator(string chrName,int start,int end){
    ti_iter_destroy(iteratorTab);
    // -1 is substracted from the start because, for some reason that is unknown to me
    // Heng Li does that in index.c in int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
    // to the get the coordinates right
    // I do this for consistency with the command line program   
    if(start>0)
	start--;
    iteratorTab=ti_query(fpTab,chrName.c_str(),start,end); 
}

bool ReadTabix::readLine(string & line){
    int length; //useless ?
    const char *buffer;

    bool toReturn=(buffer=ti_read(fpTab,iteratorTab,&length));
    if(toReturn){
	line=string(buffer);
    }

    return toReturn;
}

