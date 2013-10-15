/*
 * testReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <memory>

#include "utils.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "SimpleVCF.h"
#include "VCFreader.h"
#include "FilterVCF.h"

#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

using namespace std;

int main (int argc, char *argv[]) {


    int minGQcutoff=40;
    int minMQcutoff=30;
    int minCovcutoff=0;
    int maxCovcutoff=1000;
    SetVCFFilters  * filtersVCF;
    bool allowCloseIndelProx = false;
    bool allowRepeatMasking  = false;
    bool allowSysErr         = false;
    bool allowall            = false;
    bool allowallMQ          = false;

    double     minMapabilitycutoff =0;//map20 is wrong, we need to fix it

    const string usage=string(string(argv[0])+
			      "This program filters VCF files (prints to the stdout)\n"+
			      +" [options] <vcf file>"+"\n"+
			      "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+

			      "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
			      "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
			      "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      "\t"+"--allowall"         +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      "\t"+"--allowallMQ"       +"\t\t" +"Allow all sites but still filter on MQ        (default: "+booleanAsString(allowallMQ)+")\n");



			      //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+

			      
    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //last arg is program name
    for(int i=1;i<(argc-1);i++){ 

	if(strcmp(argv[i],"--minCov") == 0 ){
	    minCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

       
	if(strcmp(argv[i],"--maxCov") == 0 ){
	    maxCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--minGQ") == 0 ){
	    minGQcutoff=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--minMQ") == 0 ){
	    minMQcutoff=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }
	

	if(strcmp(argv[i],"--allowindel") == 0 ){
	    allowCloseIndelProx =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowrm") == 0 ){
	    allowRepeatMasking   =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowSysErr") == 0 ){
	    allowSysErr     =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowall") == 0 ){
	    allowall     =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowallMQ") == 0 ){
	    allowallMQ     =true;
	    continue;
	}

	cerr<<"Wrong option "<<argv[i]<<endl;
	return 1;
    }

    //    cerr<<"minMQcutoff "<<minMQcutoff<<endl;

    filtersVCF= new SetVCFFilters (minGQcutoff          ,
				   minMQcutoff          ,
				   minMapabilitycutoff  ,
				   !allowCloseIndelProx ,
				   !allowRepeatMasking  ,
				   !allowSysErr         ,
				   minCovcutoff      ,
				   maxCovcutoff  ,
				   allowall,
				   allowallMQ);
    VCFreader vcfr (string(argv[argc-1]),5);
	// 75060,
	// 75070,
	// 5);
	
    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();
	if(passedFilters(toprint,filtersVCF)){
	    cout<<*toprint<<endl;
	}

    }

    cerr<<"filtersVCF finished successfully "<<rejectFiltersTally()<<endl;
    delete(filtersVCF);
    return 0;
}

