/*
 * FilterVCF
 * Date: Aug-22-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "FilterVCF.h"


//! To check whether a SimpleVCF passes the some filters
/*!
 *
 * This subroutine checks if a SimpleVCF record passes the following tests:
 * 
 *  0) The variation is not an indel
 *  1) Has A,C,G or T as reference allele
 *  2) Has A,C,G,T or . as alternative allele
 *  3) Has sequence coverage (depth) between two given values
 *  4) If the Map20 field is present, if checks if it is above a given cutoff
 *  5) Has the root mean square of the mapability score above a certain cutoff
 *  6) Has GQ field above a certain cutoff
 *  7) If not flagged as close to indel
 *  8) If not flagged as  syserr (systematic error)
 *  9) If not flagged as rm (repeat masked)\
 * 10) If the genotype is unknown (GT field = ./.)
 *
 \param smvcf               The SimpleVCF object to check
 \param minCovcutoff        The minimum coverage cutoff
 \param maxCovcutoff        The maximum coverage cutoff
 \param minMapabilitycutoff    The mapability score cutoff
 \param minMQcutoff            The cutoff for the root mean square of the mapability score above a certain cutoff
 \param minGQcutoff            The cutoff for the GQ field (genotype quality)
 \return  : True if the SimpleVCF passes all the aforementioned checks
*/

bool passedFilters(SimpleVCF * smvcf,int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minGQcutoff){
    /////////////////////////////////////////////////////////////////
    //           FILTERING INDELS AND UNDEFINED REF ALLELE         //
    /////////////////////////////////////////////////////////////////


    bool notFoundNonindel=smvcf->containsIndel() ;

    if(notFoundNonindel){ rejectIndel++;  //we don't want indels		
	return false;
    }

    if(!smvcf->isResolvedSingleBasePairREF()){ rejectREFValidREF++; 
	return false; }  //REF al  
    if(!smvcf->isResolvedSingleBasePairALT()){ rejectREFValidALT++; 
	return false;  } //ALT al




    /////////////////////////////////
    //FILTERING ON GENOTYPE QUALITY /
    /////////////////////////////////
		  


    // - are in the 2.5% tails of the coverage distribution
    int coverageREF=smvcf->getDepth();

    if(coverageREF < minCovcutoff ||
       coverageREF > maxCovcutoff ){
	rejectLOWCOV_REF++;
	return false;
    }



    // - have Map20 < 1
    if(smvcf->hasInfoField("Map20")){
	if(smvcf->getInfoField<double>("Map20") < minMapabilitycutoff){
	    rejectMap20++;
	    return false;
	}
    }






    // - have MQ< 30
    double minMQ=smvcf->getInfoField<double>("MQ");
		    
    if(minMQ < minMQcutoff){
	rejectLOWMQ++;
	return false;
    }


    // - have GQ < 40
    //float minQual =smvcf->getQual();
    float minQual =smvcf->getGenotypeQual();
					
    if(minQual< minGQcutoff){
	rejectLOWQUAL++;
	return false;
    }

		    

	    
    // - are ± 5bp of InDels
    bool isCloseToIndel=(smvcf->getCloseIndel() );	    

    if(isCloseToIndel){
	rejectCloseIndels++;
	return false;
    }


    // - are flagged as SysErr
    if(smvcf->hasInfoField("SysErr")){
	rejectSysERR++;
	return false;
    }



    // - are flagged as RM (repeat masked)
    if(smvcf->hasInfoField("RM")){
	rejectRM++;
	return false;
    }





		    
		    
    //rejecting unknown genotype (GT= ./.)
    if(smvcf->isUnresolvedGT() ){ rejectREF_unknownGeno++; 
	return false; }



    return true;		   
}

string rejectFiltersTally(){

    return 
	"rejectIndel            "+stringify(rejectIndel)+"\n"+
	"rejectREFValidREF      "+stringify(rejectREFValidREF)+"\n"+
	"rejectREFValidALT      "+stringify(rejectREFValidALT)+"\n"+
	"rejectLOWCOV_REF       "+stringify(rejectLOWCOV_REF)+"\n"+
	"rejectMap20            "+stringify(rejectMap20)+"\n"+
	"rejectLOWMQ            "+stringify(rejectLOWMQ)+"\n"+
	"rejectLOWQUAL          "+stringify(rejectLOWQUAL)+"\n"+
	"rejectCloseIndels      "+stringify(rejectCloseIndels)+"\n"+
	"rejectSysERR           "+stringify(rejectSysERR)+"\n"+
	"rejectRM               "+stringify(rejectRM)+"\n"+
	"rejectREF_unknownGeno  "+stringify(rejectREF_unknownGeno)+"\n";


}
