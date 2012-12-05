/*
 * SimpleVCF
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "SimpleVCF.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer

SimpleVCF::SimpleVCF(){

}

SimpleVCF::SimpleVCF(string line){
    unresolvedGT=false;
    homozygousREF=false;
    heterozygous=false;
    homozygousALT=false;
    closeIndel=false;

    indexGenotype=-1; 
    indexGenotypeQual=-1; 
    indexDepth=-1;    
    indexPL=-1;       

    typeOfData=1;

    //cout<<"SimpleVCF "<<line<<endl;
    fields=allTokens(line,'\t');


    chrName=                     fields[0];
    position=string2uint(        fields[1]);

    id     =                     fields[2];

    ref    =                     fields[3];
    alt    =                     fields[4];

    //for sites with multiple alt bases
    altAlleles      =     allTokens(alt,',');
    allAltResolvedSingleBasePair=true;
    for(int i=0;i<altAlleles.size();i++){
	allAltResolvedSingleBasePair&=validAltBP( altAlleles[i] ); //true if ref = A,C,G,T or .
    }
    
    
    //boolean flags for insert
    isIndel=isInsert(ref) || isInsert(alt);
    //boolean flags for a single bp in ref or alt
    resolvedSingleBasePairREF=validOneBP(ref); //true if ref = A,C,G or T
    resolvedSingleBasePairALT=validAltBP(alt); //true if ref = A,C,G,T or .

    if(fields[5] == "."){
	qual=0.0;
    }else{
	qual   = destringify<float>( fields[5] );
    }
    filter =                     fields[6];


    //INFO FIELD
    infoFieldRaw =        fields[7] ;
    infoField    =        info2map( infoFieldRaw );


    //FORMAT FIELDS
    rawFormatNames  = fields[8];
    rawFormatValues = fields[9];
    formatFieldNames  = allTokens(rawFormatNames ,':');
    formatFieldValues = allTokens(rawFormatValues,':');
    
    if(formatFieldNames.size() != formatFieldValues.size()){
	cerr<<"SimpleVCF: for line "<<line<<" the format field does not have as many fields as the values"<<endl;
	exit(1);
    }

    for(int i=0;i<formatFieldNames.size();i++){
	// cout<<"formatFieldNames["<<i<<"] "<<formatFieldNames[i]<<endl;
	if(formatFieldNames[i] == "GT"){ 
	    indexGenotype     =i; 
	    formatFieldGT=                   formatFieldValues[i]; 
	    bool determinedGenotype=false;
	    //Taken from http://www.broadinstitute.org/gatk/guide/topic?name=intro
	    if(formatFieldGT == "./."){ determinedGenotype=true; unresolvedGT=true; }
	    if(formatFieldGT == "0/0"){ determinedGenotype=true; homozygousREF=true; }
	    if(formatFieldGT == "0/1"){ determinedGenotype=true; heterozygous=true; }
	    if(formatFieldGT == "1/1"){ determinedGenotype=true; homozygousALT=true; }
	    if(formatFieldGT == "1/2"){ determinedGenotype=true; heterozygousALT=true; }

	    if(!determinedGenotype){
		cerr<<"SimpleVCF: unable to determine genotype for line "<<line<<""<<endl;
		exit(1);
	    }
	    continue;
	}

	if(formatFieldNames[i] == "GQ"){ indexGenotypeQual =i; formatFieldGQ=destringify<float>(formatFieldValues[i]); continue; }
	if(formatFieldNames[i] == "DP"){ indexDepth        =i; formatFieldDP=destringify<int>  (formatFieldValues[i]); continue;}
	if(formatFieldNames[i] == "PL"){ 
	    indexPL        = i; 
	    formatFieldPL  = formatFieldValues[i];
	    vector<string> plfields = allTokens(formatFieldPL,',');
	    if(plfields.size() != 3){
		cerr<<"SimpleVCF: for line "<<line<<" the PL field does not have 3 fields"<<endl;
		exit(1);
	    }
	    formatFieldPLHomoRef =  destringify<int>(plfields[0]);
	    formatFieldPLHetero  =  destringify<int>(plfields[1]);
	    formatFieldPLHomoAlt =  destringify<int>(plfields[2]);
	    continue;
	}

	if(formatFieldNames[i] == "A"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(int j=0;j<adfield.size();j++){
		countA.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames[i] == "C"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(int j=0;j<adfield.size();j++){
		countC.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames[i] == "G"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(int j=0;j<adfield.size();j++){
		countG.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}

	if(formatFieldNames[i] == "T"){   
	    vector<string> adfield = allTokens( formatFieldValues[i] ,',');
	    for(int j=0;j<adfield.size();j++){
		countT.push_back(   destringify<int>( adfield[j]) );
	    }
	    continue;
	}
	    


    }

    // cout<<getADforA()<<endl;
    // cout<<getADforC()<<endl;
    // cout<<getADforG()<<endl;
    // cout<<getADforT()<<endl;

    // cout<<"end"<<endl;

}

SimpleVCF::~SimpleVCF(){
    // cout<<"DESTRUCTOR SimpleVCF"<<endl;
    // delete  altAlleles;
    // delete  fields;
    // delete  infoField;
    // delete formatFieldNames;
    // delete formatFieldValues;
    //exit(1);
}


string SimpleVCF::getRef() const{
    return ref;
}

string SimpleVCF::getAlt() const{
    return alt;
}

string SimpleVCF::getID() const{
    return id;
}

string SimpleVCF::getChr() const{
    return chrName;
}

unsigned int SimpleVCF::getPosition() const{
    return position;
}

string  SimpleVCF::getFilter() const{
    return filter;
}

string SimpleVCF::getInfoFieldRaw() const{
    return infoFieldRaw;
}


bool   SimpleVCF::hasInfoField(string tag) const{
    return (infoField.find(tag)  != infoField.end());
}





string  SimpleVCF::getGenotype() const{
    if(indexGenotype != -1){
	return formatFieldGT;
    }else{
	return "";
    }
}


float   SimpleVCF::getGenotypeQual() const{
    if(indexGenotypeQual != -1){
	return formatFieldGQ;
    }else{
	return -1.0;
    }
}

int     SimpleVCF::getDepth() const{
    if(indexDepth != -1){
	return formatFieldDP;
    }else{
	return -1;
    }
}

int     SimpleVCF::getDepthInfo() const{
    if(hasInfoField("DP")){
	int toReturn = getInfoField<int>("DP");
	return toReturn;
    }
    return -1;
}

string  SimpleVCF::getPL() const{
    if(indexPL != -1){
	return formatFieldPL;
    }else{
	return "";
    }
}

int     SimpleVCF::getPLHomoRef() const{
    if(indexPL != -1){
	return formatFieldPLHomoRef;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHetero() const{
    if(indexPL != -1){
	return formatFieldPLHetero;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHomoAlt() const{
    if(indexPL != -1){
	return formatFieldPLHomoAlt;
    }else{
	return -1;
    }
}



float     SimpleVCF::getQual() const{
    return qual;
}

void    SimpleVCF::setCloseIndel(bool closeIndel){
    this->closeIndel=closeIndel;
}

bool SimpleVCF::getCloseIndel() const{
    return closeIndel;
}


bool SimpleVCF::containsIndel() const{
    return isIndel;
}


bool    SimpleVCF::isResolvedSingleBasePairREF() const{
    return resolvedSingleBasePairREF;
}

bool    SimpleVCF::isResolvedSingleBasePairALT() const{
    return resolvedSingleBasePairALT;
}


bool SimpleVCF::isUnresolvedGT() const{
    return unresolvedGT;
}

bool SimpleVCF::isHomozygousREF() const{
    return homozygousREF;
}

bool SimpleVCF::isHeterozygous() const{
    return heterozygous;
}

bool SimpleVCF::isHomozygousALT() const{
    return homozygousALT;
}

bool SimpleVCF::isHeterozygousALT() const{
    return heterozygousALT;
}

char SimpleVCF::getRandomAllele() const{
    if( homozygousREF ){ return ref[0]; }
    if( homozygousALT ){ return alt[0]; }
    if( heterozygous ){ 
	//pick an allele at random
	if(randomBool())
	    return ref[0]; 
	else
	    return alt[0]; 
    }
    //error
    //should we allow heterozygous alt ?
    cerr<<"SimpleVCF: Cannot generate a random allele for "<<(*this)<<endl;
    exit(1);    
}


// ostream& operator<<(ostream& os, const SimpleVCF& smvcf){
void SimpleVCF::print(ostream& os) const{

    // os<<smvcf.chrName<<"\t"
    //   <<smvcf.position<<"\t"
    //   <<smvcf.ref<<"\t"
    //   <<smvcf.alt<<"\t"
    //   <<smvcf.qual<<"\t"
    //   <<smvcf.filter<<"\t"
    //   <<smvcf.infoFieldRaw<<"\t"
    //   <<smvcf.rawFormatNames<<"\t"
    //   <<smvcf.rawFormatValues;

    os<<chrName<<"\t"
      <<position<<"\t"
      <<ref<<"\t"
      <<alt<<"\t"
      <<qual<<"\t"
      <<filter<<"\t"
      <<infoFieldRaw<<"\t"
      <<rawFormatNames<<"\t"
      <<rawFormatValues;


    // return os;
}


bool SimpleVCF::hasAllele(char bp) const {

    if(resolvedSingleBasePairREF && resolvedSingleBasePairALT){ //only look at sites with a single bp
	if( homozygousREF ){ return (ref[0] == bp); }
	if( homozygousALT ){ return (alt[0] == bp); }
	if( heterozygous ){ 
	    return ( (ref[0] == bp) || (alt[0] == bp));
	}
    } 
    return false;

}


bool SimpleVCF::hasAtLeastOneA() const  {   
    return hasAllele('A');
}

bool SimpleVCF::hasAtLeastOneC() const  {
    return hasAllele('C');
}

bool SimpleVCF::hasAtLeastOneG() const  {
    return hasAllele('G');
}

bool SimpleVCF::hasAtLeastOneT() const  {
    return hasAllele('T');
}


bool SimpleVCF::hasAllele(int indexAlle) const  {

    if(indexAlle == 1){return hasAllele('A');}
    if(indexAlle == 2){return hasAllele('C');}
    if(indexAlle == 3){return hasAllele('G');}
    if(indexAlle == 4){return hasAllele('T');}

    cerr<<"SimpleVCF: hasAllele() request for value: "<<indexAlle<<" cannot be completed, must be between 1 and 4, exiting"<<endl;
    exit(1);

}


pair<int,int> SimpleVCF::returnLikelyAlleleCountForRefAlt(int minPLdiffind) const{

    if(unresolvedGT) //unresolved, we cannot infer anything
	return pair<int,int>(0,0);

    if ( (formatFieldPLHetero-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind) {  //high likelihood of homo ref, produce 2 alleles ref
	// refAlleles+=2;
	// altAlleles+=0;
	return pair<int,int>(2,0);
    } else{
	if ((formatFieldPLHetero-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind) {  //high likelihood of homo alt, produce 2 alleles alt
	    // refAlleles+=0;
	    // altAlleles+=2;
	    return pair<int,int>(0,2);
	} else {
	    if ((formatFieldPLHomoRef-formatFieldPLHetero) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHetero) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
		// refAlleles+=1;
		// altAlleles+=1;
		return pair<int,int>(1,1);
	    }else{
		if ((formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoAlt) <minPLdiffind ) { //high likelihood of at least one alt, produce 1 allele alt 
		    // refAlleles+=0;
		    // altAlleles+=1;
		    return pair<int,int>(0,1);
		}else{
		    if ( (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoRef) < minPLdiffind ) { // high likelihood of at least one ref, produce 1 allele ref
			// refAlleles+=1;
			// altAlleles+=0;
			return pair<int,int>(1,0);
		    }else{
			// refAlleles+=0;
			// altAlleles+=0;
			return pair<int,int>(0,0);
		    }
		}
	    }
	}
    }// end all cases

}








char SimpleVCF::getRandomAlleleUsingPL(int minPLdiffind) const{
    
    if(unresolvedGT)
	return 'X'; //unresolved

    if ( (formatFieldPLHetero-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind) {  //high likelihood of homo ref, produce 2 alleles ref
	// cout<<position<<"\t"<<"2,0"<<endl;
	return ref[0];
    } else{
	if ((formatFieldPLHetero-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind) {  //high likelihood of homo alt, produce 2 alleles alt
	    // cout<<position<<"\t"<<"0,2"<<endl;
	    return alt[0];
	} else {
	    if ((formatFieldPLHomoRef-formatFieldPLHetero) >= minPLdiffind && (formatFieldPLHomoAlt-formatFieldPLHetero) >= minPLdiffind) { //high likelihood of hetero, produce 1 allele of each
		// cout<<position<<"\t"<<"1,1"<<endl;
		if(randomBool())
		    return ref[0]; 
		else
		    return alt[0]; 

	    }else{
		if ((formatFieldPLHomoRef-formatFieldPLHomoAlt) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoAlt) <minPLdiffind ) { //high likelihood of at least one alt, produce 1 allele alt 
		    // cout<<position<<"\t"<<"0,1"<<endl;
		    return alt[0];
		}else{
		    if ( (formatFieldPLHomoAlt-formatFieldPLHomoRef) >= minPLdiffind && (formatFieldPLHetero-formatFieldPLHomoRef) < minPLdiffind ) { // high likelihood of at least one ref, produce 1 allele ref
			// cout<<position<<"\t"<<"1,0"<<endl;
			return ref[0]; 
		    }else{
			return 'X'; //unresolved
		    }
		}
	    }
	}
    }// end all cases
}

int SimpleVCF::getADforA(){
    int toReturn=0;
    for(int j=0;j<countA.size();j++){
	toReturn+=countA[j];
    }
    return toReturn;
}
 
int SimpleVCF::getADforC(){
    int toReturn=0;
    for(int j=0;j<countC.size();j++){
	toReturn+=countC[j];
    }
    return toReturn;
}

int SimpleVCF::getADforG(){
    int toReturn=0;
    for(int j=0;j<countG.size();j++){
	toReturn+=countG[j];
    }
    return toReturn;
}
 

int SimpleVCF::getADforT(){
    int toReturn=0;
    for(int j=0;j<countT.size();j++){
	toReturn+=countT[j];
    }
    return toReturn;
}

