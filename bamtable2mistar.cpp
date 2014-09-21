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


#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

using namespace std;



void setVarsEPO(ReadTabix * rtEPO,string & epoChr,unsigned int & epoCoord,bool & cpgEPO,char & allel_ref,char & allel_chimp,char & allel_anc,bool & lineLeftEPO,string & lineFromEPO){

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


    allel_ref   = fieldsEPO[2][0];//reference allele
    allel_anc   = fieldsEPO[3][0];//inferred ancestor
    allel_chimp = fieldsEPO[4][0];//chimp;


}


int main (int argc, char *argv[]) {
    // cout<<INT_MAX<<endl;
    // bool ancAllele           = false;

    const string usage=string(string(argv[0]) +" [options] <bam table file> <name sample> <EPO alignment file>"+"\n"+
			      "\nThis program convert bam table  files into mistar (prints to the stdout)\n"+
    // "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
    // "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
    // "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
    // "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
    // "\t"+"--minPL [pl]"       +"\t\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
    // "\t"+"--useanc"           +"\t\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 

    // "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
    // "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
    // "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      // "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"
    //);
		
    //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+
			      "");

			      
    if(argc < 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //last arg is program name
    for(int i=1;i<(argc-3);i++){ 
                                                                                                                                                                                     
	// if(strcmp(argv[i],"--useanc") == 0 ){
	//     ancAllele       =true;
	//     continue;
	// }


	cerr<<"Wrong option "<<argv[i]<<endl;
	exit(1);
    }

    //cerr<<"Filter used: "<<*filtersVCF<<endl;
    //VCFreader vcfr   (string(argv[argc-3]),5);
    BAMTABLEreader btr  (string(argv[argc-3]),5);
    string namePop  = string(argv[argc-2]);
    string epoFile  = string(argv[argc-1]);
    string epoFileidx = epoFile+".tbi";
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
    char allel_ref;
    char allel_anc;

    // 75060,
    // 75070,
    // 5);
    //header

    cout<<"#MISTAR"<<endl;    
    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    cout<<"#PG:"<<programLine<<endl;
    cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
    cout<<"#DATE: "<<getDateString()<<endl;
    

    cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<namePop<<endl;
    unsigned int totalRec=0;
    unsigned int writtenRec=0;

    while(btr.hasData()){
	BAMTableObj * toprint=btr.getData();
	totalRec++;
	


	// if(passedFilters(toprint,filtersVCF)){

	if(firstLine){
	    firstLine=false;
	    rtEPO = new ReadTabix( epoFile.c_str()  , 
				   epoFileidx.c_str()  , 
				   toprint->getChr(), 
				   int(toprint->getPosition()),INT_MAX ); //the destructor should be called automatically

	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

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

	    
	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

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

	int refIdx  =base2int(allel_ref);
	// int chimpIdx=base2int(allel_chimp);

	//need an allele that is A,C,G,T
	if(!isResolvedDNA(allel_ref)){
	    continue;
	}


	//determine alternative allele
	char alt='N';
	string s="ACGT";
	for(int i=0;i<4;i++){
	    if(s[i] != allel_ref){
		if(toprint->hasAllele(i+1) ){
		    alt=s[i];
		}
	    }
	}

	// cout<<*toprint<<endl;
	// cout<<lineFromEPO<<endl;
	string chimpString;
	string ancString;
		

	//unresolved ancestral allele (A,C,G,T)
	if(!isResolvedDNA(allel_chimp)){
	    chimpString="0,0:0";					
	}
	//resolved ancestral allele
	else{
	    if(allel_chimp == allel_ref){//no diff between chimp and reference
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
	    if(allel_anc == allel_ref){//no diff between ancestor and reference
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








	if( (alt!='N' && toprint->hasOnly2Alleles(refIdx,base2int(alt)))  || // has alternative
	    (alt=='N' && toprint->hasOnlyThisAlleles(refIdx)) ){//no alt in bam table

	    int altCount=0;
		if(alt=='N'){
		    altCount=0;
		}else{
		    altCount=toprint->countAllele(base2int(alt));
		}
		//cout<<refIdx<<endl;
		writtenRec++;

		cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
		    allel_ref<<","<<
		    alt<<"\t"<<
		    chimpString<<"\t"<<
		    ancString<<"\t"<<
		    toprint->countAllele(refIdx)<<","<<
		    altCount<<
		    ":"<<(toprint->isCpg()?"1":"0")<<endl;

	}else{
	    continue;//triallelic, skip
	}

	//     //cout<<"ok "<<alt<<endl;
	//     if( isResolvedDNA(allel_chimp) ){	//chimp allele exists A,C,G,T

	// 	// if(allel_chimp != allel_ref ){ //triallelic with respect to chimp,skip
	// 	//     cout<<*toprint<<endl;
	// 	//     cout<<lineFromEPO<<endl;		   
	// 	// }
		
	// 	if(alt!= 'N')
	// 	    if(allel_chimp != allel_ref && //triallelic with respect to chimp,skip
	// 	       allel_chimp != alt){
	// 		continue;
	// 	    }
	    
	    
	// 	if(alt=='N' && 
	// 	   allel_chimp != allel_ref){//no alt in bamtable and chimp is different: chimp becomes alt
	// 	    alt=allel_chimp;		  
	// 	}


	// 	cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 	    allel_ref<<","<<
	// 	    alt<<"\t";

	// 	if(allel_chimp == allel_ref ){
	// 	    cout<<"1,0:"+string(cpgEPO?"1":"0")<<"\t";
	// 	}else{ 
	// 	    if(allel_chimp == alt){
	// 		cout<<"0,1:"+string(cpgEPO?"1":"0")<<"\t";
	// 	    }else{
	// 		cerr<<"Error, wrong state between "<<(*toprint)<<"\t"<<lineFromEPO<<endl;
	// 		return 1;	
	// 	    }
	// 	}



	// 	int altCount=0;
	// 	if(alt=='N'){
	// 	    altCount=0;
	// 	}else{
	// 	    altCount=toprint->countAllele(base2int(alt));
	// 	}
	// 	cout<<toprint->countAllele(refIdx)<<","<<
	// 	    altCount<<
	// 	    ":"<<(toprint->isCpg()?"1":"0")<<endl;
	//     }else{
	// 	int altCount=0;
	// 	if(alt=='N'){
	// 	    altCount=0;
	// 	}else{
	// 	    altCount=toprint->countAllele(base2int(alt));
	// 	}
	// 	//cout<<refIdx<<endl;
	// 	cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 	    allel_ref<<","<<
	// 	    alt<<"\t"<<
	// 	    "0,0:"+string(cpgEPO?"1":"0")<<"\t"<<
	// 	    toprint->countAllele(refIdx)<<","<<
	// 	    altCount<<
	// 	    ":"<<(toprint->isCpg()?"1":"0")<<endl;
	//     }

	// }else{
	//     continue;//triallelic, skip
	// }


	// if(isResolvedDNA(allel_chimp)){	//chimp allele exists
	//     if(allel_chimp != allel_ref){ //chimp allele differs
		
	// 	if(toprint->hasOnly2Alleles(refIdx,chimpIdx)){//fine, chimp allele is alt
	// 	    cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 		allel_ref<<","<<
	// 		allel_chimp<<"\t"<<
	// 		"0,1:"+string(cpgEPO?"1":"0")<<"\t"<<
	// 		toprint->countAllele(refIdx)<<","<<
	// 		toprint->countAllele(chimpIdx)<<
	// 		":"<<(toprint->isCpg()?"1":"0")<<endl;	
	// 	}else{//triallelic
	// 	    continue;
	// 	}

	//     }
	//     //reference and chimp are equal
	//     else{

	// 	if( toprint->hasOnlyThisAlleles(refIdx) ){//no alternative allele
		    
	// 	    cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 		allel_ref<<","<<
	// 		"N"<<"\t"<<
	// 		"1,0:"+string(cpgEPO?"1":"0")<<"\t"<<
	// 		toprint->countAllele(refIdx)<<","<<
	// 		"0"<<
	// 		":"<<(toprint->isCpg()?"1":"0")<<endl;	

	// 	}else{//the BAMtable has an alternative
	// 	    char alt='N';
	// 	    string s="ACGT";
	// 	    for(int i=0;i<4;i++){
	// 		if(s[i] != allel_ref){
	// 		    if(toprint->hasAllele(i+1) ){
	// 			alt=s[i];
	// 		    }
	// 		}
	// 	    }

	// 	    if(toprint->hasOnly2Alleles(refIdx,alt)){//new alternative is alt
	// 		cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 		    allel_ref<<","<<
	// 		    alt<<"\t"<<
	// 		    "1,0:"+string(cpgEPO?"1":"0")<<"\t"<<
	// 		    toprint->countAllele(refIdx)<<","<<
	// 		    toprint->countAllele(base2int(alt))<<
	// 		    ":"<<(toprint->isCpg()?"1":"0")<<endl;				
	// 	    }else{
	// 		continue;//triallelic
	// 	    }

	// 	}
	//     }
	// }
	// //no ancestral information
	// else{
	//     //maybe BAMTABLE has alternative ?
	//     char alt='N';
	//     string s="ACGT";
	//     for(int i=0;i<4;i++){
	// 	//if bamtable has non-ref allele
	// 	if(s[i] != allel_ref){
	// 	    if(toprint->hasAllele(i+1) ){
	// 		alt=s[i];
	// 	    }
	// 	}
	//     }

	//     if(toprint->hasOnly2Alleles(refIdx,alt)){//new alternative is alt
	// 	cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 	    allel_ref<<","<<
	// 	    alt<<"\t"<<
	// 	    "0,0:"+string(cpgEPO?"1":"0")<<"\t"<<
	// 	    toprint->countAllele(refIdx)<<","<<
	// 	    toprint->countAllele(base2int(alt))<<
	// 	    ":"<<(toprint->isCpg()?"1":"0")<<endl;				
	//     }else{
	// 	continue;//triallelic
	//     }
	    
	// }


	

	

	// pair<int,int> pairCount= toprint->returnLikelyAlleleCountForRefAlt(minPLdiffind);

	// if(pairCount.first != 0 || pairCount.second != 0 ){
	//     char alt=(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
	//     string chimpString;
		
	//     // cout<<allel_chimp<<"\t"<<toprint->getRef()[0]<<endl;
		

	//     //unresolved ancestral allele (A,C,G,T)
	//     if(!isResolvedDNA(allel_chimp)){
	// 	chimpString="0,0:0";					
	//     }
	//     //resolved ancestral allele
	//     else{
	// 	if(allel_chimp == toprint->getRef()[0]){//no diff between chimp and reference
	// 	    chimpString="1,0:"+string(cpgEPO?"1":"0");
	// 	}else{
	// 	    if(alt == 'N'){//no alt defined, the chimp becomes the alt
			    
	// 		alt = allel_chimp;
	// 		chimpString="0,1:"+string(cpgEPO?"1":"0");
	// 	    }else{
	// 		if(alt == allel_chimp){//no alt defined, the chimp becomes the alt
	// 		    chimpString="0,1:"+string(cpgEPO?"1":"0");
	// 		}else{ //tri-allelic site, discard
	// 		    continue;
	// 		}
	// 	    }
	// 	}
	//     }				   
	//     cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
	// 	toprint->getRef()<<","<<
	// 	alt<<"\t"<<
	// 	chimpString<<"\t"<<
	// 	pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	
	// }
	// //<<endl;
	

    }

    //epoFileFP.close();
    //    delete(filtersVCF);
    delete(rtEPO);

    cerr<<"Program "<<argv[0]<<" looked at  "<<totalRec<<" records, wrote "<<writtenRec<<" terminated gracefully"<<endl;

    
    return 0;
}

