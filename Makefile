CXX      = g++ #-g -pg
LIBGAB   = /home/gabriel_renaud/lib/
LIBTABIX = /home/gabriel_renaud/Software/tabix-0.2.6/

CXXFLAGS = -Wall -Wunused-variable -lm -O3 -lz -I${LIBGAB} -I${LIBTABIX} -Igzstream/ -c
LDFLAGS  = -lz


all:  ReadTabix.o testReadTabix SimpleVCF.o CoreVCF.o testVCF FilterVCF.o BAMTableObj.o BAMTABLEreader.o AlleleInfoReader.o mergeBAMTable filterVCF FastQObj.o FastQParser.o testReadFastq SetVCFFilters.o vcf2mistar bamtable2mistar affyVCF2mistar 23andme2mistar testMultiReadTabix vcfMulti2mistar

#mergeBAMTable.o:	mergeBAMTable.cpp
#	${CXX} ${CXXFLAGS} mergeBAMTable.cpp

mergeBAMTable:	mergeBAMTable.o ${LIBGAB}utils.o BAMTableObj.o BAMTABLEreader.o  ReadTabix.o  ${LIBTABIX}libtabix.a gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterVCF:	filterVCF.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcf2mistar:	vcf2mistar.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

vcfMulti2mistar: 	vcfMulti2mistar.o ${LIBGAB}utils.o  MultiVCFreader.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o   BAMTableObj.o BAMTABLEreader.o FilterVCF.o SetVCFFilters.o gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)


affyVCF2mistar:	affyVCF2mistar.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a  gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

23andme2mistar:	23andme2mistar.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a  gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

bamtable2mistar:	bamtable2mistar.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a   BAMTableObj.o BAMTABLEreader.o  gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@

#ReadTabix.o:	ReadTabix.cpp
#	${CXX} ${CXXFLAGS} $^

#SetVCFFilters.o:	SetVCFFilters.cpp
#	${CXX} ${CXXFLAGS} $^

#SimpleVCF.o:	SimpleVCF.cpp
#	${CXX} ${CXXFLAGS} $^

#FastQObj.o:	FastQObj.cpp
#	${CXX} ${CXXFLAGS} $^

#FastQParser.o:	FastQParser.cpp
#	${CXX} ${CXXFLAGS} $^

#testReadFastq.o:	testReadFastq.cpp
#	${CXX} ${CXXFLAGS} $^

testReadFastq: testReadFastq.o FastQObj.o FastQParser.o ${LIBGAB}utils.o gzstream/gzstream.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

#BAMTableObj.o:	BAMTableObj.cpp
#	${CXX} ${CXXFLAGS} $^

#testReadTabix.o:	testReadTabix.cpp
#	${CXX} ${CXXFLAGS} testReadTabix.cpp

#VCFreader.o:	VCFreader.cpp
#	${CXX} ${CXXFLAGS} VCFreader.cpp

#testVCF.o:	testVCF.cpp
#	${CXX} ${CXXFLAGS} testVCF.cpp


#FilterVCF.o:	FilterVCF.cpp
#	${CXX} ${CXXFLAGS} FilterVCF.cpp

testVCF:	testVCF.o ${LIBGAB}utils.o SimpleVCF.o CoreVCF.o FilterVCF.o BAMTableObj.o SetVCFFilters.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testReadTabix:	testReadTabix.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testMultiReadTabix: testMultiReadTabix.o MultiVCFreader.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o CoreVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o  gzstream/gzstream.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f  testReadTabix mergeBAMTable filterVCF testReadFastq bamtable2mistar  vcf2mistar affyVCF2mistar 23andme2mistar bamtable2mistar testMultiReadTabix vcfMulti2mistar *.o 

