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
#include "MultiVCFreader.h"
#include "BAMTableObj.h"
#include "BAMTABLEreader.h"

using namespace std;

int main (int argc, char *argv[]) {


    
    // if(1){

    // ReadTabix rt ("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz",
    // 		  "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.10.mod.vcf.gz.tbi",
    // 		  "10",
    // 		  75060,
    // 		  75070);
    // cout<<rt.getHeader()<<endl;
    // string buffer;//=new string();
    // while(rt.readLine(buffer)){
    // 	cout<<"buffer "<<buffer<<endl;
    // 	SimpleVCF svcf (buffer);
    // 	cout<<"alt "<<svcf.getAlt()<<endl;
    // 	//cout<<svcf.getInfoField<float>("MQ")<<endl;
    // 	cout<<"pl "<<svcf.getPLHomoRef()<<endl;
    // }
    // exit(1);
    // }    


    //testing tabix
    if(0){
	MultiVCFreader vcfr ("/tmp/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
			     "/tmp/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi",
			     "5",
			     1,
			     10,
			     5);
	// 557572,
	// 557582,0);
	// 75060,
	// 75070,
	// 5);
	vcfr.repositionIterator("5",60000,60026);
	
	while(vcfr.hasData()){
	    vector<SimpleVCF *> * toprint=vcfr.getMultipleData();
	    //double testingmq= 
	    //cout<<toprint->getCloseIndel()<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getPosition()<<endl;
	    //cout<<toprint->getInfoField<double>("MQ")<<endl;
	    //cout<<toprint->getDepth()<<"\t"<<toprint->getDepthInfo()<<"\t"<<( toprint->getADforA() + toprint->getADforC() + toprint->getADforG() +toprint->getADforT() )<<"\t"<<(*toprint)<<endl;
	    //cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<"\t"<<endl;
	    //     toprint->hasAtLeastOneA()<<"\t"<<
	    //     toprint->hasAtLeastOneC()<<"\t"<<
	    //     toprint->hasAtLeastOneG()<<"\t"<<
	    //     toprint->hasAtLeastOneT()<<"\t"<<endl;	    
	}
	exit(1);
    }




    if(1){
	MultiVCFreader vcfr ("/tmp/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
			     // "/tmp/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi",
			     // "5",
			     // 1,
			     // 10,
			     5);
	// 557572,
	// 557582,0);
	// 75060,
	// 75070,
	// 5);
	// vcfr.repositionIterator("5",60000,60026);
	
	while(vcfr.hasData()){

	    vector<SimpleVCF *> * toprint=vcfr.getMultipleData();
	    //double testingmq= 
	    //cout<<toprint->getCloseIndel()<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getPosition()<<endl;
	    //cout<<toprint->getInfoField<double>("MQ")<<endl;
	    //cout<<toprint->getDepth()<<"\t"<<toprint->getDepthInfo()<<"\t"<<( toprint->getADforA() + toprint->getADforC() + toprint->getADforG() +toprint->getADforT() )<<"\t"<<(*toprint)<<endl;
	    //cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<"\t"<<endl;
	    //     toprint->hasAtLeastOneA()<<"\t"<<
	    //     toprint->hasAtLeastOneC()<<"\t"<<
	    //     toprint->hasAtLeastOneG()<<"\t"<<
	    //     toprint->hasAtLeastOneT()<<"\t"<<endl;	    
	}
	exit(1);
    }






    // if(0){



    // VCFreader vcfr2 (string(argv[1]),
    // 		     7);

    // 		    // 75060,
    // 		    // 75070,
    // 		    // 5);
    
    // while(vcfr2.hasData()){
    // 	SimpleVCF * toprint=vcfr2.getData();
    // 	cout<<*toprint<<"\t"<<toprint->containsIndel()<<"\t"<<toprint->getCloseIndel()<<endl;
    // }
    // }


    // if(0){

    // BAMTABLEreader btr ("/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz",
    // 			"/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz.tbi",
    // 			"10",
    // 			936576,
    // 			1000000);
    // 			// 4465012,
    // 			// 4465542);

    // 			//			164127,
    // 			//164240);

    

    // while(btr.hasData()){
    // 	BAMTableObj * toprint=btr.getData();
    // 	cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<endl;
    // 	// toprint->setCpg(true);
    // 	// cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<"\t"<<boolStringify(toprint->isCpg())<<"\t"<<toprint->hasAtLeastOneA()<<endl;
    // }
    // }
    // if(0){

    // AlleleInfoReader * alr = new BAMTABLEreader  ("/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz",
    // 						  "/mnt/expressions/susanna/Denisova_Molar_2_L9133/divergence_estimate/files_without_coverage_higher_than_2/L9133_merged_qc_uniq_ALL_l35q37_clipped_map1_notri_2cov.chr10.bed.gz.tbi",
    // 						  "10",
    // 						  4465012,
    // 						  4465042);

    // cout<<"#############"<<endl;
    //  while(alr->hasData()){
    // 	 AlleleInfo	  * toprint=alr->getData();
    // 	 cout<<*toprint<<"\t"<<toprint->getRandomAllele()<<endl;
    // }
    // }
    return 0;
}

