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


OrfAnalyzer::OrfAnalyzer(SignalSensors &sensors,int MIN_ORF_LEN)
  : sensors(sensors), MIN_ORF_LEN(MIN_ORF_LEN)
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



Essex::CompositeNode *OrfAnalyzer::noncodingToCoding(
				    const GffTranscript &refTrans,
				    const String &refStr,
				    const Sequence &refSeq,
				    const GffTranscript &altTrans,
				    const String &altStr,
				    const Sequence &altSeq,
				    int &refOrfLen,
				    int &altOrfLen)
{
  float refStartScore; int refGenomicStart;
  refOrfLen=altOrfLen=0;
  GffTranscript *refORF=findORF(refTrans,refStr,refSeq,refStartScore,
				refGenomicStart,refOrfLen);
  float altStartScore; int altGenomicStart;
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,altStartScore,
				altGenomicStart,altOrfLen);
  bool change=false;
  if(!refORF && altORF && altOrfLen>=MIN_ORF_LEN) change=true;
  else if(refORF && altORF && refORF<MIN_ORF_LEN && altOrfLen>=MIN_ORF_LEN &&
	  altOrfLen>=2*refOrfLen)
    change=true;
  Essex::CompositeNode *ret=change ? altORF->toEssex() : NULL;
  delete refORF; delete altORF;
  return ret;
}



Essex::CompositeNode *
OrfAnalyzer::earlierStartCodon(const GffTranscript &altTrans,
			       const String &altStr,
			       const Sequence &altSeq,
			       const String &refStr,
			       const Sequence &refSeq,
			       const CigarAlignment &alignment,
			       int &oldOrfLen,
			       int &newOrfLen,
			       float &oldStartCodonScore,
			       float &newStartCodonScore)
{
  /* Requirements: either the new start codon didn't exist in the reference,
     or it had a much weaker score, or it was in a different reading frame.
   */
  int altGenomicStart;
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,newStartCodonScore,
				altGenomicStart,newOrfLen);

}


