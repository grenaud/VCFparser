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
    string getRef();
//! Retrieves the alternative allele(s) as a string as it is in the raw VCF
    string getAlt();
//! Retrieves the alternative alleles as a vector of strings
    vector<string> getAltAlleles();
//! Retrieves the # of alternative alleles found
    int getAltCount();


//! Retrieves the potential associated ID
    string getID();
//! Retrieves the chr
    string getChr();
//! Retrieves the position on the chr
    unsigned int getPosition();
//! Retrieves the filter string 
    string getFilter();

//! Retrieves the overall quality
    float  getQual();

//! To check if one of the INFO fields exist
/*!
 *
 *
  \param tag : Tag of the field in the info line
  \return  : returns true if field exists, false otherwise
  \sa  getInfoField()
*/
    bool hasInfoField(string tag);

    char getRandomAllele();

//! To retrieve the information for one of the INFO fields 
/*!
 *
 * The correct invocation has to be used. <int> for int <float> for float etc.. 

  \param tag : Tag of the field in the info line
  \return  : The information associated with the tag, if the tag does not exist, return empty string
  \sa  hasInfoField()
*/
    template <typename T>
	T   getInfoField(string tag){

	map<string,string>::iterator it=infoField.find(tag);
	if(it == infoField.end()){
	    return T();
	}else{
	    return destringify<T>(it->second);
	}
	
    }

    //! Retrieves the info field
    string getInfoField();

    //! Retrieves the genotype
    string getGenotype();
    //! Retrieves the genotype quality
    float  getGenotypeQual();
    //! Retrieves the depth
    int    getDepth();
    //! Retrieves the pl field
    string getPL();
    //! Retrieves the homo ref likelihood from the pl field
    int    getPLHomoRef();
    //! Retrieves the heterozygous likelihood from the pl field
    int    getPLHetero();
    //! Retrieves the homo alt likelihood from the pl field
    int    getPLHomoAlt();



    //! Sets the close to indel flag
    void    setCloseIndel(bool closeIndel);
    bool    getCloseIndel();
    bool    containsIndel();

    //! true if the REF allele is A,C,G or T
    bool    isResolvedSingleBasePairREF();
    //! true if the ALT allele is A,C,G,T or .
    bool    isResolvedSingleBasePairALT();
   //! true if every ALT allele is A,C,G,T or .
    bool    areAllAltResolvedSingleBasePair();

    //! true if GT field is "./."
    bool isUnresolvedGT(); 
    //! true if GT field is "0/0"
    bool isHomozygousREF();  
    //! true if GT field is "0/1"
    bool isHeterozygous();   
    //! true if GT field is "1/1"
    bool isHomozygousALT();
   //! true if GT field is "1/2"
    bool isHeterozygousALT();

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

};
#endif
