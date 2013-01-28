/*
 * testReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <memory>

#include "utils.h"
#include "AlleleInfo.h"
#include "ReadTabix.h"
#include "SimpleVCF.h"
#include "VCFreader.h"
#include "FilterVCF.h"

#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

using namespace std;

int main (int argc, char *argv[]) {


    // VCFreader vcfr ("/mnt/scratch/sergi/exome/sidron/snps/gatk_2/annotated_LowQualDeamination/Y_sidron_exome_hg19_1000g_LowQualDeamination.md.bam.annotated_1.3-14gatk.vcf.gz",
    // 		    "/mnt/scratch/sergi/exome/sidron/snps/gatk_2/annotated_LowQualDeamination/Y_sidron_exome_hg19_1000g_LowQualDeamination.md.bam.annotated_1.3-14gatk.vcf.gz.tbi",
    // 		    "Y",
    // 		    1,
    // 		    28134309,
    //                 5);


    VCFreader vcfr (string(argv[1]),5);
	// 75060,
	// 75070,
	// 5);
	
    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();
	// if(passedFilters(toprint,0, 10000,1.0,30,0)){
	cout<<toprint->getRandomAlleleUsingPL(30)<<"\t"<<*toprint<<endl;	    
	// }else{
	//     // cout<<"FAIL\t"<<*toprint<<endl;
	// }
    	//cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<endl;
    }

    //    cout<<rejectFiltersTally()<<endl;

    return 0;
}

