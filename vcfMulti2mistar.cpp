/*
 * testReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno [at here] gmail.com
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
#include "MultiVCFreader.h"
// #include "FilterVCF.h"

#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;


void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

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



    allel_anc   = fieldsEPO[3][0];//inferred ancestor
    allel_chimp = fieldsEPO[4][0];//chimp;


}


int main (int argc, char *argv[]) {
    // cout<<INT_MAX<<endl;
    int minPLdiffind=33;
    // int minGQcutoff=0;
    // int minMQcutoff=0;

    // int minCovcutoff=0;
    // int maxCovcutoff=1000;
    // // SetVCFFilters  * filtersVCF;
    // // bool allowCloseIndelProx = false;
    // // bool allowRepeatMasking  = false;
    // // bool allowSysErr         = false;
    // // bool allowall            = false;
    // // bool allowallMQ          = false;

    // bool ancAllele           = false;

    //double     minMapabilitycutoff =0;
    string epoFile   = "none";
    bool   epoFileB  = false;

    const string usage=string(string(argv[0]) +" [options] <vcf file> "+"\n"+
			      "\nThis program convert multi VCF files into mistar (prints to the stdout)\n"+
			      // "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      // "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      // "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      // "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
			      // "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
			      "\t"+"--epo [EPO alignment file]"       +"\t\t" +"Use file as EPO alignment   (default: none)\n"+ 			      
			      "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
			      // // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 

			      // "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
			      // "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
			      // "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      // "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      // "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n"
			      ""
			      );
		
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+


    for(int i=1;i<(argc-1);i++){ 
                                                                                                                                                                                     
        if(strcmp(argv[i],"--minPL") == 0 ){
            minPLdiffind=destringify<int>(argv[i+1]);
	    //            specifiedPL  =true;
            i++;
            continue;
	}

	if(strcmp(argv[i],"--epo") == 0 ){
            epoFile=string(argv[i+1]);
	    epoFileB=true;
            i++;
            continue;
	}

	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }
			      
    if(argc < 2 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    MultiVCFreader vcfr   (string(argv[argc-1]),5);



    string epoFileidx = epoFile+".tbi";

    cerr<<"VCF file "<<(string(argv[argc-1]))<<endl;
    //cerr<<"Name pop "<<(string(argv[argc-2]))<<endl;
    if(epoFileB)
	cerr<<"EPO file "<<(epoFile)<<endl;

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
    char allel_anc;

    // 75060,
    // 75070,
    // 5);
    cout<<"#MISTAR"<<endl;    
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    

    cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<vectorToString(vcfr.getPopulationNames(),"\t")<<endl;
    //cout<<"#chr\tcoord\tREF,ALT\troot\t"<<namePop<<endl;

    while(vcfr.hasData()){
    	vector<SimpleVCF *> * toprint=vcfr.getMultipleData();
	//cerr<<"line"<<endl;
	if(toprint->at(0)->containsIndel() || //skip indels
	   (toprint->at(0)->getAltAlleleCount() > 1  )){ //skip tri-allelic
	    continue;
	}

	if(firstLine){
	    firstLine=false;

	    if(epoFileB){
		rtEPO = new ReadTabix( epoFile.c_str()  , 
				       epoFileidx.c_str()  , 
				       toprint->at(0)->getChr(), 
				       int(toprint->at(0)->getPosition()),INT_MAX ); //the destructor should be called automatically
		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	    }
	}

	if(epoFileB){

	    if(!lineLeftEPO){
		cerr<<"Error, no data in the EPO file "<< toprint->at(0)->getChr() <<":"<< int(toprint->at(0)->getPosition()) <<endl;
		return 1;
	    }

	    if(epoChr != toprint->at(0)->getChr()){
		cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<toprint->at(0)->getChr()<<endl;
		return 1;
	    }



	    while(epoCoord != toprint->at(0)->getPosition()){
		if(epoCoord > toprint->at(0)->getPosition()){
		    cerr<<"Error1, are all the sites in EPO there? Difference between coords VCF="<<(*toprint->at(0))<<"\tEPO="<<lineFromEPO<<endl;
		    return 1;
		}

		if( (toprint->at(0)->getPosition() - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO->repositionIterator(toprint->at(0)->getChr() , int(toprint->at(0)->getPosition()),INT_MAX);
		    //cout<<"repo "<<int(toprint->at(0)->getPosition())<<endl;
		}


		setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);	

	    }

	    if(epoCoord != toprint->at(0)->getPosition()){
		cerr<<"Error2, are all the sites in EPO there? Difference between coords VCF="<<( toprint->at(0)->getPosition() )<<"\tline="<<lineFromEPO<<endl;
		return 1;
	    }

	}
	
	//pair<int,int> pairCount= toprint->at(0)->returnLikelyAlleleCountForRefAlt(minPLdiffind);

	// if(pairCount.first != 0 || pairCount.second != 0 ){
	char alt=(toprint->at(0)->getAlt()=="."?'N':toprint->at(0)->getAlt()[0]);
	string chimpString;
	string ancString;
	// cout<<lineFromEPO<<endl;
	// cout<<allel_chimp<<"\t"<<allel_anc<<"\t"<<epoFileB<<"\t"<<toprint->at(0)->getRef()[0]<<"\t"<<alt<<endl;
		
	if(!epoFileB){ //no epo file
	    chimpString="0,0:0";
	    ancString="0,0:0";
	    // cout<<allel_chimp<<"\t"<<toprint->at(0)->getRef()[0]<<endl;
	}else{	

	    //unresolved ancestral allele (A,C,G,T)
	    if(!isResolvedDNA(allel_chimp)){
		chimpString="0,0:0";					
	    }
	    //resolved ancestral allele
	    else{
		if(allel_chimp == toprint->at(0)->getRef()[0]){//no diff between chimp and reference
		    chimpString="1,0:"+string(cpgEPO?"1":"0");
		}else{
		    if(alt == 'N'){//no alt defined, the chimp becomes the alt			    
			alt = allel_chimp;
			chimpString="0,1:"+string(cpgEPO?"1":"0");
		    }else{
			if(alt == allel_chimp){//alt is chimp
			    chimpString="0,1:"+string(cpgEPO?"1":"0");
			}else{ //tri-allelic site, discard
			    continue;
			}
		    }
		}
	    }


	    if(!isResolvedDNA(allel_anc)){
		ancString="0,0:0";					
	    }
	    //resolved ancestral allele
	    else{
		if(allel_anc == toprint->at(0)->getRef()[0]){//no diff between ancestor and reference
		    ancString="1,0:"+string(cpgEPO?"1":"0");
		}else{
		    if(alt == 'N'){//no alt defined, the ancestor becomes the alt			    
			alt = allel_anc;
			ancString="0,1:"+string(cpgEPO?"1":"0");			    
		    }else{
			if(alt == allel_anc){//alt is ancestor
			    ancString="0,1:"+string(cpgEPO?"1":"0");
			}else{ //tri-allelic site, discard
			    continue;
			}
		    }
		}
	    }
	}
		

	cout<<toprint->at(0)->getChr()<<"\t"<< toprint->at(0)->getPosition()<<"\t"<<
	    toprint->at(0)->getRef()<<","<<
	    alt<<"\t"<<
	    chimpString<<"\t"<<
	    ancString<<"\t";

	for(unsigned int i=0;i<toprint->size();i++){
	    pair<int,int> pairCount;
	    if(toprint->at(i)->getObservedGL() ||toprint->at(i)->getObservedPL() ){
		pairCount = toprint->at(i)->returnLikelyAlleleCountForRefAlt(minPLdiffind);
	    }else{//just use GT 
		pairCount = toprint->at(i)->returnLikelyAlleleCountForRefAltJustGT();
	    }

	    cout<<pairCount.first<<","<<pairCount.second<<":"<<(toprint->at(i)->isCpg()?"1":"0");	

	    //cout<<toprint->at(i)->getAlleCountBasedOnGT();//	pairCount.first<<","<<pairCount.second<<":"<<(toprint->at(0)->isCpg()?"1":"0")<<endl;	
	    if(i != (toprint->size()-1))
		cout<<"\t";
	}
	//	}
	cout<<endl;
	//	}

    }

    //epoFileFP.close();
    //delete(filtersVCF);
    delete(rtEPO);
    cerr<<"Program terminated gracefully"<<endl;
    return 0;
}

