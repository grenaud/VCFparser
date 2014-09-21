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
#include <gzstream.h>

#include "ReadTabix.h"
#include "utils.h"


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
    // int minPLdiffind=33;
    // int minGQcutoff=0;
    // int minMQcutoff=0;

    // int minCovcutoff=0;
    // int maxCovcutoff=1000;
    // SetVCFFilters  * filtersVCF;
    // bool allowCloseIndelProx = false;
    // bool allowRepeatMasking  = false;
    // bool allowSysErr         = false;
    // bool allowall            = false;
    // bool allowallMQ          = false;

    // // bool ancAllele           = false;

    // double     minMapabilitycutoff =0;

    const string usage=string(string(argv[0]) +" [options] <vcf file>  <EPO alignment file>"+"\n"+
			      "\nThis program convert simple affymetrix VCF files into mistar (prints to the stdout)\n");
			      // "\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			      // "\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+
			      // "\t"+"--minGQ  [gq]" +"\t\t"+"Minimal genotype quality (default: "+stringify(minGQcutoff)+")\n"+
			      // "\t"+"--minMQ  [mq]" +"\t\t"+"Minimal mapping quality (default: "+stringify(minMQcutoff)+")\n"+
			      // "\t"+"--minMap [minmap]" +"\t"+"Minimal mapability (default: "+stringify(minMapabilitycutoff)+")\n"+
	
			      // "\t"+"--minPL [pl]"       +"\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+ 
			      // // "\t"+"--useanc"           +"\t\t" +"Use inferred ancestor instead of chimp      (default: "+stringify(ancAllele)+")\n"+ 

			      // "\t"+"--allowindel"       +"\t\t" +"Allow sites considered within 5bp of an indel (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
			      // "\t"+"--allowrm"          +"\t\t" +"Allow sites labeled repeat masked             (default: "+booleanAsString(allowRepeatMasking)+")\n"+
			      // "\t"+"--allowSysErr"      +"\t\t" +"Allow sites labeled as systematic error       (default: "+booleanAsString(allowSysErr)+")\n"+
			      // "\t"+"--allowall"          +"\t\t" +"Allow all sites                               (default: "+booleanAsString(allowall)+")\n"+
			      // "\t"+"--allowallMQ"        +"\t\t" +"Allow all sites   but still filter on MQ      (default: "+booleanAsString(allowall)+")\n");
		
	      //"\t"+"--minPL  [pl]" +"\t\t"+"Use this as the minimum difference of PL values instead of GQ    (default: "+stringify(minPLdiffind)+")\n"+

			      
    if(argc != 3 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    //    cerr<<"Filter used: "<<*filtersVCF<<endl;
    string vcffile =  string(argv[argc-2]);
    //string namePop  = string(argv[argc-2]);
    string epoFile  = string(argv[argc-1]);
    string epoFileidx = epoFile+".tbi";

    cerr<<"VCF file "<<(string(argv[argc-3]))<<endl;
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
    char allel_anc;



    igzstream myVCFgzfile;
    myVCFgzfile.open(vcffile.c_str(), ios::in);

    if (!myVCFgzfile.good()){
	cerr << "Unable to open file "<<vcffile<<endl;
	return 1;
    }

    

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
    

    //cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<namePop<<endl;
    //cout<<"#chr\tcoord\tREF,ALT\troot\t"<<namePop<<endl;
    string line;
    bool hasData= getline (myVCFgzfile,line);
    vector<string> popNames;

    //header
    while (hasData ){
	if(!line.empty()){
	    // cout<<"line "<<line<<endl;
	    if(line[0] == '#'){
		if(strBeginsWith(line,"#CHROM")){
		    vector<string> vs=allTokens(line,'\t');
		    
		    //cout<<vectorToString(vs,"-")<<endl;
		    popNames = vector<string> (vs.begin() +9 ,vs.end());
		    
		    //cout<<"TEST" <<vectorToString(popNames,"*")<<endl;
		    
		}
	    }else{
		break;
	    }
	}

	hasData= getline (myVCFgzfile,line);
     }
   
    if(popNames.empty()){
	cerr<<"Cannot get pop names, must have a header to have the names of the populations"<<endl;
	return 1;
    }
    cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<vectorToString(popNames,"\t")<<endl;

    while(hasData){
	 // cout<<"line1 "<<line<<endl;
	vector<string> vsField=allTokens(line,'\t');
	string         chr = vsField[0];
	unsigned int coord = destringify<unsigned int>(vsField[1]);
	string         refS = vsField[3];
	string         altS = vsField[4];
	// cout<<"line2 "<<line<<endl;

	char alt;
	string chimpString;
	string ancString;
	vector<string> genoToPrint;

	if(altS.size() != 1)
	    goto nextline;
	
	if(firstLine){
	    firstLine=false;
	    rtEPO = new ReadTabix( epoFile.c_str()  , 
				   epoFileidx.c_str()  , 
				   chr, 
				   int(coord),INT_MAX ); //the destructor should be called automatically
	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);

	}

	// cout<<"line3 "<<line<<endl;


	if(!lineLeftEPO){
	    cerr<<"Error, no data in the EPO file "<< chr <<":"<< coord <<endl;
	    return 1;
	}

	if(epoChr != chr){
	    cerr<<"Error, the chromosome does not match the one in the EPO file = "<<epoChr <<" and not "<<chr<<endl;
	    return 1;
	}


	//	cout<<"line2 "<<line<<endl;
	while(epoCoord != coord){
	    if(epoCoord > coord){
		cerr<<"Error, are all the sites in EPO there? Difference between coords "<<line<<"\t"<<lineFromEPO<<endl;
		return 1;
	    }

	    if( (coord - epoCoord ) >= limitToReOpenFP){ //seeking in the file
		rtEPO->repositionIterator(chr , int(coord),INT_MAX);
	    }


	    setVarsEPO(rtEPO,epoChr,epoCoord,cpgEPO,allel_chimp,allel_anc,lineLeftEPO,lineFromEPO);
	}

	if(epoCoord != coord){
	    cerr<<"Error, are all the sites in EPO there? Difference between coords "<<line<<"\t"<<lineFromEPO<<endl;
	    return 1;
	}

	alt=(altS=="."?'N':altS[0]);

	// cout<<"line4 "<<line<<endl;

	//unresolved ancestral allele (A,C,G,T)
	if(!isResolvedDNA(allel_chimp)){
	    chimpString="0,0:0";					
	}
	//resolved ancestral allele
	else{
	    if(allel_chimp == refS[0]){//no diff between chimp and reference
		chimpString="1,0:"+string(cpgEPO?"1":"0");
	    }else{
		if(alt == 'N'){//no alt defined, the chimp becomes the alt			    
		    alt = allel_chimp;
		    chimpString="0,1:"+string(cpgEPO?"1":"0");
		}else{
		    if(alt == allel_chimp){//alt is chimp
			chimpString="0,1:"+string(cpgEPO?"1":"0");
		    }else{ //tri-allelic site, discard
			goto nextline;
		    }
		}
	    }
	}

	// cout<<"line5 "<<line<<endl;

	if(!isResolvedDNA(allel_anc)){
	    ancString="0,0:0";					
	}
	//resolved ancestral allele
	else{
	    if(allel_anc == refS[0]){//no diff between ancestor and reference
		ancString="1,0:"+string(cpgEPO?"1":"0");
	    }else{
		if(alt == 'N'){//no alt defined, the ancestor becomes the alt			    
		    alt = allel_anc;
		    ancString="0,1:"+string(cpgEPO?"1":"0");			    
		}else{
		    if(alt == allel_anc){//alt is ancestor
			ancString="0,1:"+string(cpgEPO?"1":"0");
		    }else{ //tri-allelic site, discard
			goto nextline;
			//continue;
		    }
		}
	    }
	}


	

		

	cout<<chr<<"\t"<< coord<<"\t"<<
	    refS<<","<<
	    alt<<"\t"<<
	    chimpString<<"\t"<<
	    ancString<<"\t";
	    //pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	



	for(unsigned int i=9;i<vsField.size();i++){

	    //homo ref
	    if(vsField[i] == "0/0"){
		genoToPrint.push_back("2,0:0");
		continue;
	    }

	    //homo alt
	    if(vsField[i] == "1/1"){
		genoToPrint.push_back("0,2:0");
		continue;
	    }

	    //hetero
	    if(vsField[i] == "0/1" || 
	       vsField[i] == "1/0"){
		genoToPrint.push_back("1,1:0");
		continue;
	    }

	    //unresolved
	    if(vsField[i] == "./."){
		genoToPrint.push_back("0,0:0");
		continue;
	    }
	    

	    cerr<<"ERROR with symbol "<<vsField[i] <<" line "<<line<<endl;
	    exit(1);

	}

	cout<<vectorToString(genoToPrint,"\t")<<endl;
	
	nextline:
	
	hasData= getline (myVCFgzfile,line);
    }

    
    // 	SimpleVCF * toprint=vcfr.getData();



    //if(passedFilters(toprint,filtersVCF)){


    	    //<<endl;
    //}

    //}

    //epoFileFP.close();
    // delete(filtersVCF);
    delete(rtEPO);
    cerr<<"Program terminated gracefully"<<endl;
    return 0;
}

