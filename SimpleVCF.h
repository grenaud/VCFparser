/*
 * SimpleVCF
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef SimpleVCF_h
#define SimpleVCF_h

#include <string> 
#include <vector> 
#include <map>

#include "AlleleInfo.h"
#include "utils.h"

using namespace std;

//! A  class to hold VCF info for a single individual
/*!
The constructor parses the line from the VCF and populates the fields. 
*/
class SimpleVCF : public AlleleInfo{
private:
  

    //Taken from http://www.broadinstitute.org/gatk/guide/topic?name=intro
    bool unresolvedGT;   //if GT == "./."
    bool homozygousREF;  //if GT == "0/0"
    bool heterozygous;   //if GT == "0/1"
    bool homozygousALT;  //if GT == "1/1"
    bool heterozygousALT;  //if GT == "1/2"

    bool resolvedSingleBasePairREF;
    bool resolvedSingleBasePairALT;
    bool allAltResolvedSingleBasePair;

    bool closeIndel;
    bool isIndel;
    //Format fields
    string rawFormatNames;
    string rawFormatValues;

    int indexGenotype; //GT
    int indexGenotypeQual; //GQ
    int indexDepth;    //DP
    int indexPL;       //PL
    
    //Format fields
    string formatFieldGT;
    float  formatFieldGQ;
    int    formatFieldDP;
    
    //todo
    string  formatFieldPL;
    int     formatFieldPLHomoRef;
    int     formatFieldPLHetero;
    int     formatFieldPLHomoAlt;



    unsigned int position;
    string    chrName;
    string    id;

    string    ref;
    string    alt;
    vector<string> altAlleles;

    float    qual;    
    string filter;

    vector<string> fields;
    string infoFieldRaw;
    map<string, string> infoField;

    vector<string> formatFieldNames;
    vector<string> formatFieldValues;

    inline bool hasAllele(char bp) const ;
public:

//! Constructor 
/*!
  \param line : The raw line from the VCF file
*/
    SimpleVCF(string line);

    //! Dummy constructor, do not use 
    SimpleVCF();

    ~SimpleVCF();
//! Retrieves the reference allele as a string as it is in the raw VCF
    string getRef() const;
//! Retrieves the alternative allele(s) as a string as it is in the raw VCF
    string getAlt() const;
//! Retrieves the alternative alleles as a vector of strings
    vector<string> getAltAlleles() const;
//! Retrieves the # of alternative alleles found
    int getAltCount() const;


//! Retrieves the potential associated ID
    string getID() const;
//! Retrieves the chr
    string getChr() const;
//! Retrieves the position on the chr
    unsigned int getPosition() const;
//! Retrieves the filter string 
    string getFilter() const;

//! Retrieves the overall quality
    float  getQual() const;

//! To check if one of the INFO fields exist
/*!
 *
 *
  \param tag : Tag of the field in the info line
  \return  : returns true if field exists, false otherwise
  \sa  getInfoField()
*/
    bool hasInfoField(string tag) const;

    char getRandomAllele() const;
    char getRandomAlleleUsingPL(int minPLdiffind) const;


//! To retrieve the information for one of the INFO fields 
/*!
 *
 * The correct invocation has to be used. <int> for int <float> for float etc.. 

  \param tag : Tag of the field in the info line
  \return  : The information associated with the tag, if the tag does not exist, return empty string
  \sa  hasInfoField()
*/
    template <typename T>
	T   getInfoField(string tag) {

	map<string,string>::iterator it=infoField.find(tag);
	if(it == infoField.end()){
	    return T();
	}else{
	    return destringify<T>(it->second);
	}
	
    }

    //! Retrieves the info field
    string getInfoField() const; 

    //! Retrieves the genotype
    string getGenotype() const;
    //! Retrieves the genotype quality
    float  getGenotypeQual() const;
    //! Retrieves the depth
    int    getDepth() const;
    //! Retrieves the pl field
    string getPL() const;
    //! Retrieves the homo ref likelihood from the pl field
    int    getPLHomoRef() const;
    //! Retrieves the heterozygous likelihood from the pl field
    int    getPLHetero() const;
    //! Retrieves the homo alt likelihood from the pl field
    int    getPLHomoAlt() const;



    //! Sets the close to indel flag
    void    setCloseIndel(bool closeIndel);
    bool    getCloseIndel() const;
    bool    containsIndel() const;

    //! true if the REF allele is A,C,G or T
    bool    isResolvedSingleBasePairREF() const;
    //! true if the ALT allele is A,C,G,T or .
    bool    isResolvedSingleBasePairALT() const;
   //! true if every ALT allele is A,C,G,T or .
    bool    areAllAltResolvedSingleBasePair() const;

    //! true if GT field is "./."
    bool isUnresolvedGT() const; 
    //! true if GT field is "0/0"
    bool isHomozygousREF() const;  
    //! true if GT field is "0/1"
    bool isHeterozygous() const;   
    //! true if GT field is "1/1"
    bool isHomozygousALT() const;
   //! true if GT field is "1/2"
    bool isHeterozygousALT() const;

    /* friend ostream& operator<<(ostream& os, const SimpleVCF& smvcf); */
    void print(ostream& os) const;

    //! true if the current record contains at least one A
    bool hasAtLeastOneA() const  ;
    //! true if the current record contains at least one C
    bool hasAtLeastOneC() const  ;
    //! true if the current record contains at least one G
    bool hasAtLeastOneG() const  ;
    //! true if the current record contains at least one T
    bool hasAtLeastOneT() const  ;
    
    
    //returns true if this record has this allele (A:1,C:2,G:3,T:4)
    bool hasAllele(int indexAlle) const;

    //! returns the count of reference and alternative allele based on PL field
    pair<int,int> returnLikelyAlleleCountForRefAlt(int minPLdiffind=50) const;

};
#endif