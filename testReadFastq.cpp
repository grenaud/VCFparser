#include <iostream>
#include <fstream>
#include <memory>

#include "utils.h"
#include "FastQParser.h"
#include "FastQObj.h"

using namespace std;

int main (int argc, char *argv[]) {
    

    FastQParser fqp (argv[1]);

    while(fqp.hasData()){

	cout<<"l1"<<endl;
    	FastQObj * test	=fqp.getData();
    	cout<<*(test->getSeq())<<endl;
    	cout<<*(test->getQual())<<endl;
	cout<<"l2"<<endl;	
    }

    return 0;
}

