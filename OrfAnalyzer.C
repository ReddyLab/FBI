/****************************************************************
 OrfAnalyzer.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "OrfAnalyzer.H"
#include "BOOM/DnaAlphabet.H"
#include "BOOM/CodonIterator.H"
#include "BOOM/Constants.H"
#include "SignalPrinter.H"
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
  transcript->forgetCDS();
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
  Sequence seq(transcript,DnaAlphabet::global());
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
				    int &altOrfLen,
				    bool reverseStrand,
				    int altSeqLen,
				    float &refStartScore,
				    float &altStartScore,
				    String &msg)
{
  refStartScore=float(NEGATIVE_INFINITY);
  altStartScore=float(NEGATIVE_INFINITY);
  refOrfLen=altOrfLen=0;
  int refGenomicStart, altGenomicStart;
  GffTranscript *refORF=findORF(refTrans,refStr,refSeq,refStartScore,
				refGenomicStart,refOrfLen);
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,altStartScore,
				altGenomicStart,altOrfLen);
  bool change=false;
  if(!refORF && altORF && altOrfLen>=MIN_ORF_LEN)
    { change=true; msg="ref-no-start-codon"; }
  else if(refORF && altORF && refORF<MIN_ORF_LEN && altOrfLen>=MIN_ORF_LEN &&
	  altOrfLen>=2*refOrfLen)
    { change=true; msg="ref-ORF-too-short"; }
  else if(refORF && altORF && altOrfLen>=MIN_ORF_LEN && altOrfLen>=refOrfLen)
    { change=true; msg="possible-misannotation"; }
  if(change) {
    altORF->computePhases();
    altORF->loadSequence(altStr);
    if(reverseStrand) altORF->reverseComplement(altSeqLen);
  }
  Essex::CompositeNode *ret=change ? altORF->toEssex() : NULL;
  delete refORF; delete altORF;
  return ret;
}



GffTranscript *
OrfAnalyzer::earlierStartCodon(const GffTranscript &refTrans,
			       const String &refStr,
			       const Sequence &refSeq,
			       const GffTranscript &altTrans,
			       const String &altStr,
			       const Sequence &altSeq,
			       const CigarAlignment &altToRef,
			       int &oldOrfLen,
			       int &newOrfLen,
			       float &oldStartCodonScore,
			       float &newStartCodonScore,
			       String &oldStartStr,
			       String &newStartStr,
			       bool reverseStrand,
			       int altSeqLen,
			       Essex::CompositeNode *&msg)
{
  /* Requirements: either the new start codon didn't exist in the reference,
     or it had a much weaker score, or it was in a different reading frame.
   */

  msg=NULL;

  // See if there is an upstream start codon
  int altGenomicStart;
  GffTranscript *altORF=findORF(altTrans,altStr,altSeq,newStartCodonScore,
				altGenomicStart,newOrfLen);
  if(!altORF) return NULL;
  oldOrfLen=altTrans.getCDSlength();
  int oldBegin, oldEnd;
  altTrans.getCDSbeginEnd(oldBegin,oldEnd);
  if(altGenomicStart>=oldBegin) { delete altORF; return NULL; }
  SignalSensor *sensor=sensors.startCodonSensor;

  // If that upstream start codon was previously intronic, it's a good bet
  // that this is a functional change
  const int refGenomicStart=altToRef.mapApproximate(altGenomicStart,DIR_RIGHT);
  const int refLocalStart=refTrans.mapToTranscriptCoords(refGenomicStart);
  String reason;
  bool change=refLocalStart<0; // -1 means unmapped due to being intronic
  if(refLocalStart<0) reason="was-intronic";

  // If this start codon existed but had a poor score in the reference,
  // it may indeed indicate a functional change
  const int offset=sensor->getConsensusOffset();
  const int windowLen=sensor->getContextWindowLength();
  if(refLocalStart>=0) { // it was exonic in the ref annotation
    GffTranscript refCopy(refTrans);
    refCopy.loadSequence(refStr);
    String refRNA=refCopy.getFullSequence();
    const int begin=refLocalStart-offset;
    if(!sensor->consensusOccursAt(refRNA,refLocalStart)) {
      change=true;
      reason="bad-consensus";
    }
    else if(begin>=0) {
      Sequence rnaSeq(refRNA,DnaAlphabet::global());
      float refScore=sensor->getLogP(rnaSeq,refRNA,begin);
      if(refScore<sensor->getCutoff()) {
	change=true;
	reason="score-below-threshold";
      }
    }
  }

  // If no change, just return
  if(!change) return NULL;

  // Compute oldStartCodonScore 
  const int altLocal=altTrans.mapToTranscriptCoords(oldBegin);
  GffTranscript altCopy(altTrans);
  altCopy.loadSequence(altStr);
  String altRNA=altCopy.getFullSequence();
  Sequence altRNAseq(altRNA,DnaAlphabet::global());
  oldStartCodonScore=sensor->getLogP(altRNAseq,altRNA,altLocal-offset);
  oldStartStr=SignalPrinter::print(*sensor,altLocal-offset,altRNA);
  const int newAltLocal=altTrans.mapToTranscriptCoords(altGenomicStart);
  newStartStr=SignalPrinter::print(*sensor,newAltLocal-offset,altRNA);

  // If the start codon existed previously but was in a different frame,
  // it's again worth reporting as possibly impacting function
  /*int refBegin, refEnd;
  refTrans.getCDSbeginEnd(refBegin,refEnd);
  const int localRefBegin=refTrans.mapToTranscriptCoords(refBegin);
  int oldFrame=(localRefBegin-refLocalStart)%3;
  const int oldLocal=altTrans.mapToTranscriptCoords(oldBegin);
  const int newLocal=altTrans.mapToTranscriptCoords(altGenomicStart);
  int newFrame=(oldLocal-newLocal)%3;
  if(newFrame==0 && oldFrame!=0) {
    change=true;
    cout<<"case D"<<endl;
  }
  //if(newFrame!=0 && oldFrame==0) change=true;
  */

  msg=new Essex::CompositeNode("reference");
  msg->append(reason);
  altORF->computePhases();
  altORF->loadSequence(altStr);
  return altORF;
}


