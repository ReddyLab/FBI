/****************************************************************
 OrfAnalyzer.H
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_OrfAnalyzer_H
#define INCL_OrfAnalyzer_H
#include <iostream>
#include "BOOM/GffTranscript.H"
#include "BOOM/String.H"
#include "BOOM/Sequence.H"
#include "SignalSensors.H"
using namespace std;
using namespace BOOM;


/*  This class works only on forward-strand genes.
 */
class OrfAnalyzer {
public:
  OrfAnalyzer(SignalSensors &);
  GffTranscript *findORF(const GffTranscript &,const String &genomeStr,
			 const Sequence &genomeSeq,float &startCodonScore,
			 int &genomicStartCodonPos,int &orfLength);
private:
  SignalSensors &sensors;
  int findStartCodon(const String &transcript,float &startCodonScore);
};

#endif

