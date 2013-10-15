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
#include <gzstream.h>

#include "utils.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"


// #include "BAMTableObj.h"
// #include "BAMTABLEreader.h"

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

    const string usage=string(string(argv[0]) +" [options] <23andme file> <name sample> <EPO alignment file>"+"\n"+
			      "\nThis program convert 23andme files into mistar (prints to the stdout)\n"+
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
    // BAMTABLEreader btr  (),5);
    string inputFile23 = string(argv[argc-3]);
    igzstream myfile;
    myfile.open(inputFile23.c_str(), ios::in);

    if (!myfile.good()){
        cerr << "Unable to open file "<<inputFile23<<endl;
        return 1;
    }

    
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
    string line;
    bool hasData= getline (myfile,line);
    while (hasData ){

	// while(btr.hasData()){
	// BAMTableObj * toprint=btr.getData();
	vector<string> token= allTokens(line,'\t');
	// cout<<"line "<<line<<endl;
	string chr    = token[1];
	unsigned int pos = destringify<unsigned int>(token[2]);
	string genotype    = token[3];

	// if(passedFilters(toprint,filtersVCF)){

	if(firstLine){
	    firstLine=false;
	    rtEPO = new ReadTabix( epoFile.c_str()  , 
				   epoFileidx.c_str()  , 
				   chr,
				   int(pos),INT_MAX ); //the destructor should be called automatically

	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_ref,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

	}


	if(!lineLeftEPO){
	    cerr<<"Error, no data in the EPO file "<< chr <<":"<< int(pos) <<endl;
	    return 1;
	}

	if(epoChr != chr){
	    cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<chr<<endl;
	    return 1;
	}


	// cout<<"2"<<endl;
	while(epoCoord != pos){
	    if(epoCoord > pos){
		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<(line)<<"\t"<<lineFromEPO<<endl;
		return 1;
	    }

	    if( (pos - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		rtEPO->repositionIterator(chr , int(pos),INT_MAX);
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


	if(epoCoord != pos){
	    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<line<<"\t"<<lineFromEPO<<endl;
		return 1;
	}

	//int refIdx  =base2int(allel_ref);
	// int chimpIdx=base2int(allel_chimp);
	char alt='N';
	string s="ACGT";
	string chimpString;
	string ancString;

	//need an allele that is A,C,G,T
	//cout<<pos<<"\t"<<allel_ref<<endl;
	if(!isResolvedDNA(allel_ref)){
	    goto nextline;
	}

	// cout<<"3"<<endl;

	//determine alternative allele
	for(int i=0;i<4;i++){
	    if(s[i] != allel_ref){
		if(genotype.find(s[i]) != string::npos){ //new non-ref allele found
		    //toprint->hasAllele(i+1) ){
		    alt=s[i];
		}
	    }
	}

	// cout<<*toprint<<endl;
	// cout<<lineFromEPO<<endl;
		

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
			//continue;
			goto nextline;
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
			//continue;
			goto nextline;
		    }
		}
	    }
	}


	// cout<<"4"<<endl;





	if( (alt!='N' && genotype.find(alt) != string::npos) )   // has alternative
	    {//no alt in bam table

	    int refCount=0;
	    int altCount=0;
	    string rr = stringify(allel_ref)+ stringify(allel_ref);
	    string ra = stringify(allel_ref)+ stringify(alt);
	    string ar = stringify(alt)      + stringify(allel_ref);
	    string aa = stringify(alt)      + stringify(alt);

	    if(genotype == rr ){
	    	refCount=2;
	    }else{
	    	if(genotype == aa){
	    	    altCount=2;
	    	}else{
	    	    if( (genotype == ra) ||
	    		(genotype == ar) ){
	    		refCount=1;
	    		altCount=1;
	    	    }else{
	    		//goto nextline;
			cout<<"Error: Site not found "<<line<<"\t"<<lineFromEPO<<endl;

	    	    }			
	    	}
	    }

	    //cout<<refIdx<<endl;
	    cout<<chr<<"\t"<< pos<<"\t"<<
		allel_ref<<","<<
		alt<<"\t"<<
		chimpString<<"\t"<<
		ancString<<"\t"<<
		refCount<<","<<
		altCount<<
		":"<<("0")<<endl;//no cpg info

	}else{
	    goto nextline;
	    //continue;//triallelic, skip
	}
	
        nextline:
        
        hasData = getline (myfile,line);
	
    }

    //epoFileFP.close();
    //    delete(filtersVCF);
    delete(rtEPO);
    return 0;
}

