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

    ReadTabix rt ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
    		  "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
    		  "10",
    		  75060,
    		  75070);
    string buffer;//=new string();
    while(rt.readLine(buffer)){
    	cout<<"buffer "<<buffer<<endl;
    	SimpleVCF svcf (buffer);
    	cout<<svcf.getAlt()<<endl;
    	//cout<<svcf.getInfoField<float>("MQ")<<endl;
    	cout<<svcf.getPLHomoRef()<<endl;
    }
    }    
    
    if(1){
    VCFreader vcfr ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
    		    "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
    		    "10",
    		    557572,
    		    557582,
    		    5);

    		    // 75060,
    		    // 75070,
    		    // 5);
    
    while(vcfr.hasData()){
    	SimpleVCF * toprint=vcfr.getData();
	//double testingmq= 
	cout<<toprint->getInfoField<double>("MQ")<<endl;
	//cout<<toprint->getInfoField<double>("MQ")<<endl;
	//cout<<toprint->getDepth()<<"\t"<<toprint->getDepthInfo()<<"\t"<<( toprint->getADforA() + toprint->getADforC() + toprint->getADforG() +toprint->getADforT() )<<"\t"<<(*toprint)<<endl;
	//cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<"\t"<<endl;
	//     toprint->hasAtLeastOneA()<<"\t"<<
	//     toprint->hasAtLeastOneC()<<"\t"<<
	//     toprint->hasAtLeastOneG()<<"\t"<<
	//     toprint->hasAtLeastOneT()<<"\t"<<endl;

    }
    }




    if(0){



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


    if(0){

    BAMTABLEreader btr ("/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz",
    			"/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz.tbi",
    			"10",
			936576,
			1000000);
			// 4465012,
			// 4465542);

    			//			164127,
    			//164240);

    

    while(btr.hasData()){
    	BAMTableObj * toprint=btr.getData();
    	cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<endl;
	// toprint->setCpg(true);
    	// cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<"\t"<<toprint->hasAtLeastOneA()<<endl;
    }
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
	 cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<endl;
    }
    }
    return 0;
}

