#include <iostream>
#include <fstream>
#include <memory>

#include "utils.h"
#include "FastQParser.h"
#include "FastQObj.h"

using namespace std;

int main (int argc, char *argv[]) {
    

    //FastQParser fqp (argv[1],true);
    FastQParser fqp (argv[1]);

    while(fqp.hasData()){
    	FastQObj * test	=fqp.getData();
	cout<<*(test)<<endl;
    	// cout<<*(test->getSeq())<<endl;
    	// cout<<*(test->getQual())<<endl;
    }



    return 0;
}

