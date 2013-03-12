/*
 * testVCF
 * Date: Aug-21-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "utils.h"
#include "AlleleInfo.h"
#include "SimpleVCF.h"
#include "FilterVCF.h"
#include "SetVCFFilters.h"

#include "BAMTableObj.h"

using namespace std;

int main (int argc, char *argv[]) {
    

    SimpleVCF * test =new SimpleVCF ("20	1242215	.	A	T	892.16	.	AC=1,1;AF=0.50,0.50;AN=2;BaseQRankSum=1.855;DP=78;Dels=0.00;FS=1.011;HRun=7;HaplotypeScore=37.9652;MQ=38.64;MQ0=0;MQRankSum=-1.420;QD=11.44;ReadPosRankSum=2.138;Map20=1.0,1000gALT=G;AF1000g=0.04;AMR_AF=0.04;EUR_AF=0.09;TS=HPGOMC;TSseq=A,-,-,-,-,-;CAnc=-;GAnc=-;OAnc=-;bSC=897	GT:DP:GQ:PL:A:C:G:T:IR	1/1:78:40:922,0,1586:0,2:0,0:12,34:10,20:1");    


   SetVCFFilters filtersToUse (40,//gq
			       20, //mq
			       1.0, //mapability
			       true, //filter indel
			       true, //RM
			       true, //syserr
			       70,    //cov min
			       77, //cov max
			       false,
			       false);  //let all sites through

   cout<<boolStringify(passedFilters(test,  &filtersToUse))<<endl;


   //QA performed Gabriel Jan 21 2013

   //GQ ok,      passed
   //GQ too low, passed

   //MQ ok,      passed
   //MQ too low, passed

   //Map20 fail,      passed
   //Map20 passed,    passed

   //indel was check using debugging

   //rm fail,         passed
   //rm passed,       passed
   //rm do not check, passed

   //syserr fail,         passed
   //syserr passed,       passed
   //syserr do not check, passed


   //cov min, ok          passed
   //cov min, too low     passed


   //cov max, ok          passed
   //cov max, too low     passed


    return 0;

    SimpleVCF * test4 =new SimpleVCF ("20	1242215	.	A	G,T	892.16	.	AC=1,1;AF=0.50,0.50;AN=2;BaseQRankSum=1.855;DP=78;Dels=0.00;FS=1.011;HRun=7;HaplotypeScore=37.9652;MQ=38.64;MQ0=0;MQRankSum=-1.420;QD=11.44;ReadPosRankSum=2.138;1000gALT=G;AF1000g=0.04;AMR_AF=0.04;EUR_AF=0.09;RM;TS=HPGOMC;TSseq=A,-,-,-,-,-;CAnc=-;GAnc=-;OAnc=-;bSC=897	GT:DP:GQ:PL:A:C:G:T:IR	1/2:78:99:922,0,1586:0,2:0,0:12,34:10,20:1");    

    SimpleVCF * test2 =new SimpleVCF ("20	1242215	.	A	T	892.16	.	AC=1,1;AF=0.50,0.50;AN=2;BaseQRankSum=1.855;DP=78;Dels=0.00;FS=1.011;HRun=7;HaplotypeScore=37.9652;MQ=38.64;MQ0=0;MQRankSum=-1.420;QD=11.44;ReadPosRankSum=2.138;1000gALT=G;AF1000g=0.04;AMR_AF=0.04;EUR_AF=0.09;TS=HPGOMC;TSseq=A,-,-,-,-,-;CAnc=-;GAnc=-;OAnc=-;bSC=897	GT:DP:GQ:PL:A:C:G:T:IR	1/1:78:99:922,0,1586:0,2:0,0:12,34:10,20:1");

    SimpleVCF * test3 =new SimpleVCF("10	464212	.	A	C	161.39	.	AC=1;AF=0.50;AN=2;BaseQRankSum=0.433;DP=18;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=34.11;MQ0=1;MQRankSum=-1.010;QD=8.97;ReadPosRankSum=-1.395;RM;TS=HGO;TSseq=A,-,G;GAnc=A;OAnc=G;mSC=0.002;pSC=0.004	GT:DP:GQ:PL:A:C:G:T:IR	0/1:18:99:191,0,175:7,1:8,1:0,0:0,0:0");


    BAMTableObj bo1 ("10	137643	137644	0	1	0	0");
    BAMTableObj bo2 ("10	94313	94314	0	1	0	0");

    cout<<bo1<<endl;
    cout<<bo2<<endl;
    cout<<(bo1+bo1)<<endl;

    cout<<test->isHomozygousREF()<<endl;
    cout<<test->isHeterozygous()<<endl;
    cout<<test->isHomozygousALT()<<endl;
    cout<<test->isHeterozygousALT()<<endl;

    //bool passedFilters(SimpleVCF smvcf,int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minQCcutoff);
    SetVCFFilters filtersToUse2 (20,//gq
				20, //mq
				1.0,				
				true,
				true,
				true,
				0,
				100);
	//                              minCovcutoffREF, maxCovcutoffREF, minMapabilitycutoff, minMQcutoff, minGQcutoff);
    // 	cout<<boolStringify(passedFilters(&test,  0,     100,             1.0,                 20,            20))<<endl;
    // cout<<boolStringify(passedFilters(&test2, 4 ,     98,              1.0,                 20,           100))<<endl;
    cout<<boolStringify(passedFilters(test,  &filtersToUse2))<<endl;
    // cout<<boolStringify(passedFilters(&test2, 4 ,     98,              1.0,                 20,           100))<<endl;


    cout<<rejectFiltersTally()<<endl;

    cout<<"random = "<<test3->getRandomAllele()<<endl;
    cout<<"random = "<<test3->getRandomAllele()<<endl;
    cout<<"random = "<<test3->getRandomAllele()<<endl;
    cout<<"random = "<<test3->getRandomAllele()<<endl;

    vector<AlleleInfo *> vecAll;
    vecAll.push_back(test2);
    vecAll.push_back(test3);
    vecAll.push_back(&bo1);
    vecAll.push_back(&bo2);

    for(int i=0;i<vecAll.size();i++){
	cout<<vecAll[i]->getChr()<<"\t"<<vecAll[i]->getPosition()<<"\t"<<"\t"<<vecAll[i]->getRandomAllele()<<endl;
	cout<<*(vecAll[i])<<endl;
    }



    return 0;
}

