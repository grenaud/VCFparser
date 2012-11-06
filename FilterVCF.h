/*
 * FilterVCF
 * Date: Aug-22-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef FilterVCF_h
#define FilterVCF_h

#include "SimpleVCF.h"

static int rejectIndel=0;
static int rejectREFValidREF=0;
static int rejectREFValidALT=0;
static int rejectLOWCOV_REF=0;
static int rejectMap20=0;
static int rejectLOWMQ=0;
static int rejectLOWQUAL=0;
static int rejectCloseIndels=0;
static int rejectSysERR=0;
static int rejectRM=0;
static int rejectREF_unknownGeno=0;

bool passedFilters(SimpleVCF * smvcf,int minCovcutoff,int maxCovcutoff,double minMapabilitycutoff,int minMQcutoff,int minGQcutoff);
string rejectFiltersTally();

#endif
