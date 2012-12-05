CXX      = g++ -g -pg
LIBGAB   = /home/gabriel_renaud/lib/
LIBTABIX = /home/gabriel_renaud/Software/tabix-0.2.6/

CXXFLAGS = -lm -O3 -lz -I${LIBGAB} -I${LIBTABIX}  -c
LDFLAGS  = -lz


all:  ReadTabix.o testReadTabix SimpleVCF.o testVCF FilterVCF.o BAMTableObj.o BAMTABLEreader.o AlleleInfoReader.o mergeBAMTable filterVCF

mergeBAMTable.o:	mergeBAMTable.cpp
	${CXX} ${CXXFLAGS} mergeBAMTable.cpp

mergeBAMTable:	mergeBAMTable.o ${LIBGAB}utils.o BAMTableObj.o BAMTABLEreader.o  ReadTabix.o  ${LIBTABIX}libtabix.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterVCF:	filterVCF.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o FilterVCF.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

ReadTabix.o:	ReadTabix.cpp
	${CXX} ${CXXFLAGS} $^

SimpleVCF.o:	SimpleVCF.cpp ${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} $^

BAMTableObj.o:	BAMTableObj.cpp ${LIBGAB}utils.o 
	${CXX} ${CXXFLAGS} $^

testReadTabix.o:	testReadTabix.cpp
	${CXX} ${CXXFLAGS} testReadTabix.cpp

VCFreader.o:	VCFreader.cpp
	${CXX} ${CXXFLAGS} VCFreader.cpp

testVCF.o:	testVCF.cpp
	${CXX} ${CXXFLAGS} testVCF.cpp


FilterVCF.o:	FilterVCF.cpp
	${CXX} ${CXXFLAGS} FilterVCF.cpp

testVCF:	testVCF.o ${LIBGAB}utils.o SimpleVCF.o FilterVCF.o BAMTableObj.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testReadTabix:	testReadTabix.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o 
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f  testReadTabix ReadTabix.o SimpleVCF.o testReadTabix.o VCFreader.o FilterVCF.o BAMTableObj.o BAMTABLEreader.o AlleleInfoReader.o mergeBAMTable.o mergeBAMTable filterVCF.o filterVCF

