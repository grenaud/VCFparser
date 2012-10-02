CXX      = g++ -g #-pg
LIBGAB   = /home/gabriel_renaud/lib/
LIBTABIX = /home/gabriel_renaud/Software/tabix-0.2.6/

CXXFLAGS = -lm -O3 -lz -I${LIBGAB} -I${LIBTABIX}  -c
LDFLAGS  = -lz


all:  ReadTabix.o testReadTabix SimpleVCF.o testVCF FilterVCF.o BAMTableObj.o BAMTABLEreader.o AlleleInfoReader.o


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
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

testReadTabix:	testReadTabix.o ${LIBGAB}utils.o  ReadTabix.o  ${LIBTABIX}libtabix.a SimpleVCF.o VCFreader.o  BAMTableObj.o BAMTABLEreader.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean :
	rm -f  testReadTabix ReadTabix.o SimpleVCF.o testReadTabix.o VCFreader.o FilterVCF.o BAMTableObj.o BAMTABLEreader.o AlleleInfoReader.o

