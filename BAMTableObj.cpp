/*
 * BAMTableObj
 * Date: Aug-28-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "BAMTableObj.h"

BAMTableObj::BAMTableObj(){
    chrName="";
    position=0;
    totalAllelCount=0;

    for(int i = 0;i<4;i++){
	alleleCount.push_back( 0 );
    }

}


BAMTableObj::BAMTableObj(string line){
    fields=allTokens(line,'\t');
    typeOfData=1;

    if(fields.size() != 7){
	cerr<<"BAMTableObj: line "<<line<<" does not have 7 fields"<<endl;
	exit(1);
    }

    timeval time;
    gettimeofday(&time, NULL);

    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );

    chrName=                     fields[0];
    position=string2uint(        fields[2]); //the first column is the 0-base coordinate
    if( (position - string2uint(fields[1])) != 1){
	cerr<<"BAMTableObj: for line "<<line<<" the first position column is not 1 minus the second"<<endl;
	exit(1);
    }

    totalAllelCount=0;    

    for(int i = 0;i<4;i++){
	int toadd=destringify<int>( fields[3+i] );
	if(toadd < 0){
	    cerr<<"BAMTableObj: for line "<<line<<" cannot have negative allele counts"<<endl;
	    exit(1);
	}
	totalAllelCount+=toadd;
	alleleCount.push_back( toadd );
    }


    if(totalAllelCount == 0){
	cerr<<"BAMTableObj: for line "<<line<<" all allele counts are 0"<<endl;
	exit(1);
    }

    
}

BAMTableObj::~BAMTableObj(){

}


unsigned int BAMTableObj::getPosition() const{
    return position;
}

string BAMTableObj::getChr() const{
    return chrName;
}

//to implement
char BAMTableObj::getRandomAllele() const{
    int randIndex=rand()%totalAllelCount;//returns a number between 0 and (totalAllelCount-1)
    for(int i = 0;i<4;i++){
	randIndex-=alleleCount[i];// if we find the allele count in the 
	if(randIndex<0){
	    switch(i){
	    case 0:
		return 'A';
	    case 1:
		return 'C';
	    case 2:
		return 'G';
	    case 3:
		return 'T';
	    default:
		cerr<<"BAMTableObj: Cannot generate a random allele for "<<(*this)<<endl;
		exit(1);    		
	    }
	}
    }

    cerr<<"BAMTableObj: Cannot generate a random allele for "<<(*this)<<endl;
    exit(1);    
}





//ostream& operator<<(ostream& os, const BAMTableObj& bto){
void BAMTableObj::print(ostream& os) const{

    // os<<bto.chrName<<"\t"
    //   <<(bto.position-1)<<"\t"
    //   <<bto.position<<"\t"
    //   <<bto.alleleCount[0]<<"\t"
    //   <<bto.alleleCount[1]<<"\t"
    //   <<bto.alleleCount[2]<<"\t"
    //   <<bto.alleleCount[3];
    
    os<<chrName<<"\t"
      <<(position-1)<<"\t"
      <<position<<"\t"
      <<alleleCount[0]<<"\t"
      <<alleleCount[1]<<"\t"
      <<alleleCount[2]<<"\t"
      <<alleleCount[3];
    

    // return os;
}

bool BAMTableObj::hasAtLeastOneA() const  {
    return (alleleCount[0] >0);
}

bool BAMTableObj::hasAtLeastOneC() const  {
    return (alleleCount[1] >0);
}

bool BAMTableObj::hasAtLeastOneG() const  {
    return (alleleCount[2] >0);
}

bool BAMTableObj::hasAtLeastOneT() const  {
    return (alleleCount[3] >0);
}


bool BAMTableObj::hasAllele(int indexAlle) const{

    if( (indexAlle<1)  ||  (indexAlle>4) ) {	
	cerr<<"BAMTableObj: hasAllele() index cannot be "<<indexAlle<<" for this record  "<<(*this)<<endl;
	exit(1);    
    }
    return ( alleleCount[indexAlle-1] > 0 );
}


int BAMTableObj::countAllele(int indexAlle) const{

    if( (indexAlle<1)  ||  (indexAlle>4) ) {	
	cerr<<"BAMTableObj: countAllele() index cannot be "<<indexAlle<<" for this record  "<<(*this)<<endl;
	exit(1);    
    }
    return ( alleleCount[indexAlle-1]  );
}




bool BAMTableObj::hasOnly2Alleles(int firstIndex,int secondIndex) const{    

    if( (firstIndex<1)  ||  (firstIndex>4) ) {	
	cerr<<"BAMTableObj: hasOnly2Alleles() index cannot be "<<firstIndex<<" for this record  "<<(*this)<<endl;
	exit(1);    
    }

    if( (secondIndex<1)  ||  (secondIndex>4) ) {	
	cerr<<"BAMTableObj: hasOnly2Alleles() index cannot be "<<secondIndex<<" for this record  "<<(*this)<<endl;
	exit(1);    
    }


    for(int i=0;i<4;i++){	
    //checking for alleles excluding the firstIndex secondIndex
	if( (i==(firstIndex-1)) ||
	    (i==(secondIndex-1))  )
	    continue;
    
	if( alleleCount[i] > 0 )
	    return false;
    }
    return true;
}



bool BAMTableObj::hasOnlyThisAlleles(int firstIndex) const{    

    if( (firstIndex<1)  ||  (firstIndex>4) ) {	
	cerr<<"BAMTableObj: hasOnlyThisAlleles index cannot be "<<firstIndex<<" for this record  "<<(*this)<<endl;
	exit(1);    
    }


    for(int i=0;i<4;i++){	
	//checking for alleles excluding the firstIndex secondIndex
	if( (i==(firstIndex-1))  )
	    continue;
    
	if( alleleCount[i] > 0 )
	    return false;
    }

    return true;
}




