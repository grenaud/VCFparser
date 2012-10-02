/*
 * SimpleVCF
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "SimpleVCF.h"


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

    // cout<<"SimpleVCF "<<line<<endl;
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
	}

	if(formatFieldNames[i] == "GQ"){ indexGenotypeQual =i; formatFieldGQ=destringify<float>(formatFieldValues[i]);}
	if(formatFieldNames[i] == "DP"){ indexDepth        =i; formatFieldDP=destringify<int>  (formatFieldValues[i]);}
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
	}
    }

    

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


string SimpleVCF::getRef(){
    return ref;
}

string SimpleVCF::getAlt(){
    return alt;
}

string SimpleVCF::getID(){
    return id;
}

string SimpleVCF::getChr(){
    return chrName;
}

unsigned int SimpleVCF::getPosition(){
    return position;
}

string  SimpleVCF::getFilter(){
    return filter;
}

string SimpleVCF::getInfoField(){
    return infoFieldRaw;
}


bool   SimpleVCF::hasInfoField(string tag){
    return (infoField.find(tag)  != infoField.end());
}




string  SimpleVCF::getGenotype(){
    if(indexGenotype != -1){
	return formatFieldGT;
    }else{
	return "";
    }
}


float   SimpleVCF::getGenotypeQual(){
    if(indexGenotypeQual != -1){
	return formatFieldGQ;
    }else{
	return -1.0;
    }
}

int     SimpleVCF::getDepth(){
    if(indexDepth != -1){
	return formatFieldDP;
    }else{
	return -1;
    }
}

string  SimpleVCF::getPL(){
    if(indexPL != -1){
	return formatFieldPL;
    }else{
	return "";
    }
}

int     SimpleVCF::getPLHomoRef(){
    if(indexPL != -1){
	return formatFieldPLHomoRef;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHetero(){
    if(indexPL != -1){
	return formatFieldPLHetero;
    }else{
	return -1;
    }
}

int     SimpleVCF::getPLHomoAlt(){
    if(indexPL != -1){
	return formatFieldPLHomoAlt;
    }else{
	return -1;
    }
}



float     SimpleVCF::getQual(){
    return qual;
}

void    SimpleVCF::setCloseIndel(bool closeIndel){
    this->closeIndel=closeIndel;
}

bool SimpleVCF::getCloseIndel(){
    return closeIndel;
}


bool SimpleVCF::containsIndel(){
    return isIndel;
}


bool    SimpleVCF::isResolvedSingleBasePairREF(){
    return resolvedSingleBasePairREF;
}

bool    SimpleVCF::isResolvedSingleBasePairALT(){
    return resolvedSingleBasePairALT;
}


bool SimpleVCF::isUnresolvedGT(){
    return unresolvedGT;
}

bool SimpleVCF::isHomozygousREF(){
    return homozygousREF;
}

bool SimpleVCF::isHeterozygous(){
    return heterozygous;
}

bool SimpleVCF::isHomozygousALT(){
    return homozygousALT;
}

bool SimpleVCF::isHeterozygousALT(){
    return heterozygousALT;
}

//to implement
char SimpleVCF::getRandomAllele(){
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
