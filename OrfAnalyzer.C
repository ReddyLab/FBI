/****************************************************************
 OrfAnalyzer.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "OrfAnalyzer.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/CodonIterator.H"
using namespace std;
using namespace BOOM;


OrfAnalyzer::OrfAnalyzer(SignalSensors &sensors)
  : sensors(sensors)
{
}



GffTranscript *OrfAnalyzer::findORF(const GffTranscript &original,
				    const String &genomeStr,
				    const Sequence &genomeSeq,
				    float &startCodonScore,
				    int &genomicStartPos,
				    int &orfLength)
{
  GffTranscript *transcript=new GffTranscript(original);
  transcript->loadSequence(genomeStr);
  String RNA=original.getFullSequence();
  int splicedStartPos=findStartCodon(RNA,startCodonScore);
  if(splicedStartPos<0) { delete transcript; return NULL; }
  genomicStartPos=transcript->mapToGenomicCoords(splicedStartPos);
  transcript->splitUTRandCDS(genomeStr,genomicStartPos,sensors.stopCodons);
  orfLength=transcript->getCDSlength();
  return transcript;
}



int OrfAnalyzer::findStartCodon(const String &transcript,
				float &startCodonScore)
{
  SignalSensor *sensor=sensors.startCodonSensor;
  const int footprint=sensor->getContextWindowLength();
  const int offset=sensor->getConsensusOffset();
  const int L=transcript.length();
  const float cutoff=sensor->getCutoff();
  Sequence seq(transcript,PureDnaAlphabet::global());
  const int last=L-footprint;
  for(int pos=0 ; pos<last ; ++pos) {
    if(!sensor->consensusOccursAt(transcript,pos+offset)) continue;
    startCodonScore=sensor->getLogP(seq,transcript,pos);
    if(startCodonScore>=cutoff) return pos+offset;
  }
  return -1;
}




