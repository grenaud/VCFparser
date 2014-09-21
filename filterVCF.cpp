/*
 * testReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno [at here] gmail.com
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


    int minGQcutoff=0;
    int minMQcutoff=0;
    int minCovcutoff=0;
    int maxCovcutoff=1000;
    SetVCFFilters  * filtersVCF;

    bool filterCloseIndelProx = false;
    bool filterRepeatMasking  = false;
    bool filterSysErr         = false;
    // bool allowall            = false;
    // bool allowallMQ          = false;

    double     minMapabilitycutoff =0;//map20 is wrong, we need to fix it

    const string usage=string(string(argv[0])+
			      "This program filters VCF files (prints to the stdout)\n"+
			      +" [options] <vcf file>"+"\n"+
			      "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+

			      "\t"+"--filterindel"       +"\t\t" +"Filter sites considered within 5bp of an indel (default: "+booleanAsString(filterCloseIndelProx)+")\n"+
			      "\t"+"--filterrm"          +"\t\t" +"Filter sites labeled repeat masked             (default: "+booleanAsString(filterRepeatMasking)+")\n"+
			      "\t"+"--filterSysErr"      +"\t\t" +"Filter sites labeled as systematic error       (default: "+booleanAsString(filterSysErr)+")\n"+
			      // "\t"+"--allowall"         +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      // "\t"+"--allowallMQ"       +"\t\t" +"Allow all sites but still filter on MQ        (default: "+booleanAsString(allowallMQ)+")\n"
			      "");



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
	

	if(strcmp(argv[i],"--filterindel") == 0 ){
	    filterCloseIndelProx =true;
	    continue;
	}

	if(strcmp(argv[i],"--filterrm") == 0 ){
	    filterRepeatMasking   =true;
	    continue;
	}

	if(strcmp(argv[i],"--filterSysErr") == 0 ){
	    filterSysErr     =true;
	    continue;
	}

	// if(strcmp(argv[i],"--allowall") == 0 ){
	//     allowall     =true;
	//     continue;
	// }

	// if(strcmp(argv[i],"--allowallMQ") == 0 ){
	//     allowallMQ     =true;
	//     continue;
	// }

	cerr<<"Wrong option "<<argv[i]<<endl;
	return 1;
    }

    //    cerr<<"minMQcutoff "<<minMQcutoff<<endl;

    filtersVCF= new SetVCFFilters (minGQcutoff          ,
				   minMQcutoff          ,
				   minMapabilitycutoff  ,
				   filterCloseIndelProx ,
				   filterRepeatMasking  ,
				   filterSysErr         ,
				   minCovcutoff      ,
				   maxCovcutoff  ,
				   false,
				   false);
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

