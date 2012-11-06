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
#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

using namespace std;

int main (int argc, char *argv[]) {


    if(0){

    ReadTabix rt ("/mnt/454/Altaiensis/users/gabriel/EPOindices/chr13.epo.gz",
    		  "/mnt/454/Altaiensis/users/gabriel/EPOindices/chr13.epo.gz.tbi",
    		  "13",
		  1,
		  30001);
    rt.repositionIterator("13",20000,30001);
    string buffer;//=new string();
    bool lineLeftEPO = rt.readLine(buffer);


    while(lineLeftEPO){
    	cout<<"buffer "<<buffer<<endl;
	lineLeftEPO = rt.readLine(buffer);
    }



    // return 0;
    ReadTabix rt2 ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
    		  "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
    		  "10",
    		  75060,
    		  75070);
    string buffer2;//=new string();
    while(rt2.readLine(buffer2)){
    	cout<<"buffer2 "<<buffer2<<endl;
    	SimpleVCF svcf (buffer2);
    	cout<<svcf.getAlt()<<endl;
    	cout<<svcf.getInfoField<float>("MQ")<<endl;
    	cout<<svcf.getPLHomoRef()<<endl;
    }
    
    
    VCFreader vcfr ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
    		    "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
    		    "10",
    		    557572,
    		    557592,
    		    5);

    		    // 75060,
    		    // 75070,
    		    // 5);
    
    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();
    	cout<<*toprint<<"\t"<<//toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<"\t"<<
	    toprint->hasAtLeastOneA()<<"\t"<<
	    toprint->hasAtLeastOneC()<<"\t"<<
	    toprint->hasAtLeastOneG()<<"\t"<<
	    toprint->hasAtLeastOneT()<<"\t"<<endl;

    }


    VCFreader vcfr2 (string(argv[1]),
    		     7);

    		    // 75060,
    		    // 75070,
    		    // 5);
    
    while(vcfr2.hasData()){
    	SimpleVCF * toprint=vcfr2.getData();
    	cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<endl;
    }
    

    }
    
    BAMTABLEreader btr ("/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz",
    			"/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz.tbi",
    			"10",
			936576,
			1000000);
			// 4465012,
			// 4465542);

    			//			164127,
    			//164240);

    btr.repositionIterator("10",936600,996845);

    while(btr.hasData()){
    	BAMTableObj * toprint=btr.getData();
    	cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<endl;
	// toprint->setCpg(true);
    	// cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<"\t"<<toprint->hasAtLeastOneA()<<endl;
    }



    
    if(0){

    AlleleInfoReader * alr = new BAMTABLEreader  ("/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz",
						  "/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz.tbi",
						  "10",
						  4465012,
						  4465042);

    cout<<"#############"<<endl;
     while(alr->hasData()){
	 AlleleInfo	  * toprint=alr->getData();
	 cout<<*toprint<<endl;
    }

    }
    return 0;
}

