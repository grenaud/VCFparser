/*
 * BAMTableObj
 * Date: Aug-28-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef BAMTableObj_h
#define BAMTableObj_h

#include <string>
#include <vector>
#include <sys/time.h> //for srand

#include "AlleleInfo.h"
#include "utils.h"

using namespace std;

//! A  class to hold the data from BAMTable
/*!
The constructor parses the line from the file and populates the fields. 
*/
class BAMTableObj : public AlleleInfo{
private:
    vector<string> fields;
    unsigned int position;
    string    chrName;
    vector<int> alleleCount; // A=index 0, C=index 1, G=index 2, T=index 3
    int totalAllelCount;
public:

//! Constructor 
/*!
  \param line : The raw line from the BAMTable file
*/
    BAMTableObj();
    BAMTableObj(string line);
    /* BAMTableObj(string   &chrName,unsigned int position,vector<int> & alleleCount,int totalAllelCount);     */
    ~BAMTableObj();

//! Retrieves the chr
    string getChr() const;
//! Retrieves the position on the chr
    unsigned int getPosition() const;

    char getRandomAllele() const;

    /* friend ostream& operator<<(ostream& os, const BAMTableObj& bto); */
    void print(ostream& os) const;


    bool hasAtLeastOneA() const  ;
    bool hasAtLeastOneC() const  ;
    bool hasAtLeastOneG() const  ;
    bool hasAtLeastOneT() const  ;

    //returns true if this record has this allele (A:1,C:2,G:3,T:4)
    bool hasAllele(int indexAlle) const;
    //returns allele count of this allele (A:1,C:2,G:3,T:4)
    int countAllele(int indexAlle) const;

    //! Determines if a BAMTable record only has no other alleles than these two indices (A:1,C:2,G:3,T:4)
    bool hasOnly2Alleles(int firstIndex,int secondIndex) const;
    //! Determines if a BAMTable record only has no other alleles than this index (A:1,C:2,G:3,T:4)
    bool hasOnlyThisAlleles(int firstIndex) const;

    friend BAMTableObj operator+(const BAMTableObj &bo1, const BAMTableObj &bo2){
	BAMTableObj sum (bo1);	
	
	if( ( bo1.getPosition() != bo2.getPosition() ) ||
	    ( bo1.getChr()      != bo2.getChr() ) ){
	    cerr<<"BAMTableObj: operator+ both operands do not have the same chr and coord, exiting"<<endl;
	    exit(1);
	}

	for(int i = 0;i<4;i++){
	    sum.alleleCount[i]  += bo2.alleleCount[i];
	    sum.totalAllelCount += bo2.alleleCount[i];
	}
	return sum;
    }



    friend void operator+=(BAMTableObj &bo1,const BAMTableObj &bo2){
	//first operand is empty
	if( ( bo1.position     == 0 ) &&
	    ( bo1.chrName      == "" ) ){
	    bo1.position = bo2.position;
	    bo1.chrName  = bo2.chrName;
	}else{

	    if( ( bo1.getPosition() != bo2.getPosition() ) ||
		( bo1.getChr()      != bo2.getChr() ) ){
		cerr<<"BAMTableObj: operator+= both operands do not have the same chr and coord, exiting"<<endl;
		exit(1);
	    }

	}

	for(int i = 0;i<4;i++){
	    bo1.alleleCount[i]  += bo2.alleleCount[i];
	    bo1.totalAllelCount += bo2.alleleCount[i];
	}
    }


};





#endif
