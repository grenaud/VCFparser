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
    BAMTableObj(string line);
    ~BAMTableObj();

//! Retrieves the chr
    string getChr();
//! Retrieves the position on the chr
    unsigned int getPosition();

    char getRandomAllele();

    /* friend ostream& operator<<(ostream& os, const BAMTableObj& bto); */
    void print(ostream& os) const;

    bool hasAtLeastOneA() const  ;
    bool hasAtLeastOneC() const  ;
    bool hasAtLeastOneG() const  ;
    bool hasAtLeastOneT() const  ;



};
#endif
