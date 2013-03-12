/*
 * testReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <climits>

#include "utils.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "SimpleVCF.h"
#include "VCFreader.h"
#include "FilterVCF.h"

#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;


void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp, bool & lineLeftEPO,string & lineFromEPO,const bool ancAllele){

    lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
    if(!lineLeftEPO){
	cerr<<"Error, missing data in the EPO file"<<endl;
	exit(1);
    }

    vector<string> fieldsEPO  = allTokens(lineFromEPO,'\t');
    epoChr                   = fieldsEPO[0];
    epoCoord                 = string2uint(fieldsEPO[1]);					
    if(fieldsEPO[9] == "1")
	cpgEPO=true;		    
    else
	cpgEPO=false;		    
    if(ancAllele){
	allel_chimp = fieldsEPO[3][0];//inferred ancestor
    }else{
	allel_chimp = fieldsEPO[4][0];//chimp;
    }

}


int main (int argc, char *argv[]) {
    // cout<<INT_MAX<<endl;
    int minPLdiffind=33;
    int minGQcutoff=0;
    int minMQcutoff=0;

    int minCovcutoff=0;
    int maxCovcutoff=1000;
    SetVCFFilters  * filtersVCF;
    bool allowCloseIndelProx = false;
    bool allowRepeatMasking  = false;
    bool allowSysErr         = false;
    bool allowall            = false;
    bool allowallMQ          = false;

    bool ancAllele           = false;

    double     minMapabilitycutoff =0;

    const string usage=string(string(argv[0]) +" [options] <vcf file> <name sample> <EPO alignment file>"+"\n"+
			      "\nThis program convert VCF files into mistar (prints to the stdout)\n"+			     
			      "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
			      "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
	
			      "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
			      "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 

			      "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
			      "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
			      "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n");
		
	      //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+

			      
    if(argc < 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //last arg is program name
    for(int i=1;i<(argc-3);i++){ 
                                                                                                                                                                                     
        if(strcmp(argv[i],"--minPL") == 0 ){
            minPLdiffind=destringify<int>(argv[i+1]);
	    //            specifiedPL  =true;
            i++;
            continue;
	}

        if(strcmp(argv[i],"--minMap") == 0 ){
            minMapabilitycutoff=destringify<double>(argv[i+1]);
            i++;
            continue;
	}
                 
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
	    allowall       =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowallMQ") == 0 ){
	    allowallMQ      =true;
	    continue;
	}

	if(strcmp(argv[i],"--useanc") == 0 ){
	    ancAllele       =true;
	    continue;
	}

	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }

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

    cerr<<"Filter used: "<<*filtersVCF<<endl;
    VCFreader vcfr   (string(argv[argc-3]),5);
    string namePop  = string(argv[argc-2]);
    string epoFile  = string(argv[argc-1]);
    string epoFileidx = epoFile+".tbi";

    cerr<<"VCF file "<<(string(argv[argc-3]),5)<<endl;
    cerr<<"Name pop "<<(string(argv[argc-2]))<<endl;
    cerr<<"EPO file "<<(string(argv[argc-1]))<<endl;

    // string epoFileidx  = epoFileidx+".tbi";
    // if(!isFile(epoFileidx)){
    // 	cerr<<"Tabix file epoFileidx not found"<<endl;
    // 	return 1;
    // }
    // igzstream epoFileFP;
    // epoFileFP.open(epoFile.c_str(), ios::in);    // open the streams
    // string epoLine;
    string epoChr;
    unsigned int epoCoord;
    //char epoREF;

    // if (epoFileFP.good()) {
    // 	//fine
    // }else{
    // 	cerr<<"Unable to open the file "<<epoFile<<endl;
    // 	exit(1);
    // }


    // getline(epoFileFP,epoLine);
    //read first line
    // if(1){
    // 	vector<string> fieldsEPO= allTokens(epoLine,'\t');
    // 	epoChr=fieldsEPO[0];
    // 	epoCoord=string2uint(fieldsEPO[1]);	
    // }
    ReadTabix * rtEPO ;
    string lineFromEPO;
    bool lineLeftEPO;
    bool cpgEPO=false;
    bool firstLine=true;
    char allel_chimp;

    // 75060,
    // 75070,
    // 5);
    cout<<"#chr\tcoord\tREF,ALT\troot\t"<<namePop<<endl;

    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();



	if(passedFilters(toprint,filtersVCF)){

	    if(firstLine){
		firstLine=false;
		rtEPO = new ReadTabix( epoFile.c_str()  , 
				       epoFileidx.c_str()  , 
				       toprint->getChr(), 
				       int(toprint->getPosition()),INT_MAX ); //the destructor should be called automatically
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,lineLeftEPO,lineFromEPO,ancAllele);

		// lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
		// vector<string> fieldsEPO = allTokens(lineFromEPO,'\t');
		// epoChr                   = fieldsEPO[0];
		// epoCoord                 = string2uint(fieldsEPO[1]);	
                // if(fieldsEPO[9] == "1")
		//     cpgEPO=true;
		// else
		//     cpgEPO=false;
		// if(ancAllele){
		//     allel_chimp = fieldsEPO[3][0];//inferred ancestor
		// }else{
		//     allel_chimp = fieldsEPO[4][0];//chimp;
		// }		    
	    }


	    if(!lineLeftEPO){
		cerr<<"Error, no data in the EPO file "<< toprint->getChr() <<":"<< int(toprint->getPosition()) <<endl;
		return 1;
	    }

	    if(epoChr != toprint->getChr()){
		cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->getChr()<<endl;
		return 1;
	    }



	    while(epoCoord != toprint->getPosition()){
		if(epoCoord > toprint->getPosition()){
		    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(*toprint)<<"\t"<<lineFromEPO<<endl;
		    return 1;
		}

		if( (toprint->getPosition() - epoCoord ) >= limitToReOpenFP){ //seeking in the file
                    rtEPO->repositionIterator(toprint->getChr() , int(toprint->getPosition()),INT_MAX);
                }

		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,lineLeftEPO,lineFromEPO,ancAllele);
		// lineLeftEPO=(rtEPO->readLine( lineFromEPO ));
		// vector<string> fieldsEPO = allTokens(lineFromEPO,'\t');
		// epoChr                   = fieldsEPO[0];
		// epoCoord                 = string2uint(fieldsEPO[1]);					
		// if(fieldsEPO[9] == "1")
		//     cpgEPO=true;		    
		// else
		//     cpgEPO=false;		    
		// if(ancAllele){
		//     allel_chimp = fieldsEPO[3][0];//inferred ancestor
		// }else{
		//     allel_chimp = fieldsEPO[4][0];//chimp;
		// }

		// if(!lineLeftEPO){
		//     cerr<<"Error, missing data in the EPO file"<<*toprint<<endl;
		//     return 1;
		// }
	    }

	    if(epoCoord != toprint->getPosition()){
		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(*toprint)<<"\t"<<lineFromEPO<<endl;
		return 1;
	    }
	    pair<int,int> pairCount= toprint->returnLikelyAlleleCountForRefAlt(minPLdiffind);

	    if(pairCount.first != 0 || pairCount.second != 0 ){
		char alt=(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
		string chimpString;
		
		// cout<<allel_chimp<<"\t"<<toprint->getRef()[0]<<endl;
		

		//unresolved ancestral allele (A,C,G,T)
		if(!isResolvedDNA(allel_chimp)){
		    chimpString="0,0:0";					
		}
		//resolved ancestral allele
		else{
		    if(allel_chimp == toprint->getRef()[0]){//no diff between chimp and reference
			chimpString="1,0:"+string(cpgEPO?"1":"0");
		    }else{
			if(alt == 'N'){//no alt defined, the chimp becomes the alt
			    
			    alt = allel_chimp;
			    chimpString="0,1:"+string(cpgEPO?"1":"0");
			}else{
			    if(alt == allel_chimp){//no alt defined, the chimp becomes the alt
				chimpString="0,1:"+string(cpgEPO?"1":"0");
			    }else{ //tri-allelic site, discard
				continue;
			    }
			}
		    }
		}				   
		cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
		    toprint->getRef()<<","<<
		    alt<<"\t"<<
		    chimpString<<"\t"<<
		    pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	
	    }
	    //<<endl;
	}

    }

    //epoFileFP.close();
    delete(filtersVCF);
    delete(rtEPO);
    return 0;
}

