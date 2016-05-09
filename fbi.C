/****************************************************************
 fbi.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
#include "BOOM/FastaReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/CombinationIterator.H"
#include "BOOM/Essex.H"
#include "BOOM/Regex.H"
#include "BOOM/ConfigFile.H"
#include "BOOM/AminoAlphabet.H"
#include "BOOM/BandedSmithWaterman.H"
#include "Labeling.H"
#include "ProjectionChecker.H"
#include "EnumerateAltStructures.H"
#include "SignalSensors.H"
#include "GarbageCollector.H"
#include "StartCodonFinder.H"
#include "NMD.H"
#include "Variant.H"
#include "VariantClassifier.H"
#include "OrfAnalyzer.H"
using namespace std;
using namespace BOOM;


static const char *PROGRAM_NAME="find-variant-signals";
static const char *VERSION="1.0";
Alphabet alphabet;

class FBI {
public:
  FBI();
  ~FBI();
  int main(int argc,char *argv[]);
private:
  NMD nmd;
  FastaWriter fastaWriter;
  SignalSensors sensors;
  String labelingFile;
  String substrate, altDefline, xmlFilename, globalCoord;
  int MAX_SPLICE_SHIFT, MIN_EXON_LEN, MIN_INTRON_LEN, NMD_DISTANCE_PARM;
  bool allowExonSkipping, allowIntronRetention, allowCrypticSites;
  bool reverseStrand, quiet;
  String CIGAR;
  OrfAnalyzer *orfAnalyzer;
  Essex::CompositeNode *root;
  Essex::Node *startCodonMsg;
  Vector<Variant> variants;
  Regex warningsRegex, errorsRegex, variantRegex, coordRegex;
  int VCFwarnings, VCFerrors;
  int refSeqLen, altSeqLen;
  float openPenalty, extendPenalty;
  int bandwidth;
  CigarAlignment *alignment, *revAlignment;
  SubstitutionMatrix<float> *substMatrix; // protein matrix
  GarbageIgnorer garbageCollector;
  GffTranscript *loadGff(const String &filename);
  void parseNoncanonicals(const String &,Set<String> &);
  String loadSeq(const String &filename);
  String loadSeq(const String &filename,String &cigar);
  void computeLabeling(GffTranscript &,Labeling &);
  void mapLabeling(Labeling &from,Labeling &to,const CigarString &);
  void mapTranscript(GffTranscript &,const CigarString &,
		     const String &outfile,const String &genome,
		     const Sequence &genomeSeq);
  bool mapExon(GffExon &,CigarAlignment &);
  void writeProtein(const String &defline,const String &protein,
		    const String &filename);
  void append(Essex::CompositeNode *,const String &tag,const String &message);
  void append(Essex::CompositeNode *,const char *tag,const char *message);
  void append(Essex::CompositeNode *,const char *tag,const String &message);
  void append(Essex::CompositeNode *,const char *tag,int);
  void appendBrokenSignals(const TranscriptSignals *,
			   Essex::CompositeNode *status);
  void writeXML();
  void processConfig(const String &filename);
  void parseConsensusList(const String &tag,ConfigFile &,Set<String> &into);
  SignalSensor *loadModel(const String &label,ConfigFile &);
  float alignProteins(const String &refStr,const String &altStr,int &matches);
  void percentMatch(int matches,int refLen,int altLen,
		    Essex::CompositeNode *parent);
  void parseVariants(const String &,Vector<Variant> &,int substrateLen);
  Essex::CompositeNode *makeEssexVariants();
  void handleProteinFate(const AlternativeStructure &,Essex::CompositeNode *,
			 const String refProtein,const String &altSeqStr,
			 const GffTranscript *refTrans);
  int enumerateAlts(Essex::CompositeNode *altTransEssex,
		    TranscriptSignals *signals,
		    Essex::CompositeNode *status,
		    const String &altSeqStr,
		    int altSeqLen,
		    const String &refProtein,
		    const GffTranscript *refTrans,
		    GffTranscript *altTrans,
		    const String &xmlFilename,
		    Essex::CompositeNode *root,
		    ostream &osFBI);
};



int main(int argc,char *argv[])
  {
    try
      {
	FBI app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



FBI::FBI()
  : warningsRegex("/warnings=(\\d+)"), errorsRegex("/errors=(\\d+)"), 
    VCFwarnings(0), VCFerrors(0), startCodonMsg(NULL), substMatrix(NULL),
    variantRegex("(\\S+):(\\S+):(\\d+):(\\d+):([^:]*):([^:]*)"),
    coordRegex("/coord=(\\S+)"), orfAnalyzer(NULL),
    alignment(NULL), revAlignment(NULL)
  {
    // ctor
  }



FBI::~FBI()
{
  delete substMatrix;
  delete orfAnalyzer;
  delete alignment;
  delete revAlignment;
  cout<<"FBI terminated successfully"<<endl;
}



int FBI::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"ce:l:qx:");
  if(cmd.numArgs()!=6)
    throw String("\n\
fbi <fbi.config> <ref.gff> <ref.fasta> <alt.fasta> <out.gff> <out.essex>\n\
     -c = sequence has been reversed, but cigar string has not\n\
     -e N = abort if vcf errors >N\n\
     -l <file> = emit a per-nucleotide labeling for the alt sequence\n\
     -q = quiet (only report transcripts with mapping issues)\n\
     -x <file> = also emit xml\n\
  alt.fasta must have a cigar string: >ID ... /cigar=1045M3I10M7D4023M ...\n\
\n");
  const String configFile=cmd.arg(0);
  const String refGffFile=cmd.arg(1);
  const String refFasta=cmd.arg(2);
  const String altFasta=cmd.arg(3);
  const String outGff=cmd.arg(4);
  const String outFBI=cmd.arg(5);
  alphabet=DnaAlphabet::global();
  if(cmd.option('l')) labelingFile=cmd.optParm('l');
  if(cmd.option('x')) xmlFilename=cmd.optParm('x');
  reverseStrand=cmd.option('c');
  quiet=cmd.option('q');

  // Read some data from files
  processConfig(configFile);
  String refSeqStr=loadSeq(refFasta), altSeqStr=loadSeq(altFasta,CIGAR);
  const Sequence refSeq(refSeqStr,alphabet), altSeq(altSeqStr,alphabet);
  refSeqLen=refSeqStr.length(), altSeqLen=altSeqStr.length();
  GffTranscript *refTrans=loadGff(refGffFile);
  refTrans->loadSequence(refSeqStr);
  String refProtein=refTrans->getProtein();

  // Set up to generate structured output in Essex/XML
  String transcriptID=refTrans->getTranscriptId();
  String geneID=refTrans->getGeneId();
  root=new Essex::CompositeNode("report");
  append(root,"substrate",substrate);
  if(!globalCoord.empty()) append(root,"global-coords",globalCoord);
  append(root,"transcript-ID",transcriptID);
  append(root,"gene-ID",geneID);
  append(root,"vcf-warnings",VCFwarnings);
  append(root,"vcf-errors",VCFerrors);
  append(root,"alignment",CIGAR);
  //append(root,"defline",altDefline);
  Essex::CompositeNode *essexVariants=makeEssexVariants();
  root->append(essexVariants);
  refTrans->computePhases();
  if(reverseStrand) refTrans->reverseComplement(refSeqLen);
  Essex::CompositeNode *refTransEssex=refTrans->toEssex();
  refTransEssex->getTag()="reference-transcript";
  VariantClassifier classifier(variants,VariantClassifier::REF,*refTrans);
  refTransEssex->append(classifier.makeVariantsNode());
  root->append(refTransEssex);
  Essex::CompositeNode *status=new Essex::CompositeNode("status");
  root->appendChild(status);
  ofstream osFBI(outFBI.c_str());
  if(cmd.option('e') && VCFerrors>cmd.optParm('e').asInt()) {
    status->append("too-many-vcf-errors");
    if(!quiet) {
      if(!xmlFilename.empty()) writeXML();
      osFBI<<*root<<endl;
      osFBI<<"#===========================================================\n";
    }
    return -1;
  }
  
  // Check that the reference gene is well-formed
  bool noStart, noStop, PTC, badSpliceSite, referenceIsOK=true;
  String msg;
  if(!ProjectionChecker::geneIsWellFormed(*refTrans,refSeqStr,
					  noStart,noStop,PTC,badSpliceSite,
					  status,sensors,NMD_DISTANCE_PARM)) {
    if(quiet) return -1;
    status->prepend("bad-annotation");
    referenceIsOK=false;
  }

  // Compute the reference labeling
  Labeling refLab(refSeqLen);
  computeLabeling(*refTrans,refLab);

  // Make CIGAR alignment
  CigarString cigar(CIGAR);
  if(reverseStrand) cigar.reverse();
  alignment=cigar.getAlignment();
  revAlignment=alignment->invert(altSeqLen);

  // Project the reference GFF over to an alternate GFF
  mapTranscript(*refTrans,cigar,outGff,altSeqStr,altSeq);

  // Generate labeling file
  Labeling projectedLab(altSeqLen);
  mapLabeling(refLab,projectedLab,cigar);
  if(!labelingFile.empty()) {
    ofstream os(labelingFile.c_str());
    os<<projectedLab;
    os.close(); }

  // Check the projection to see if the gene might be broken
  bool mapped=false;
  if(referenceIsOK) {
    GffTranscript *altTrans=loadGff(outGff);
    altTrans->loadSequence(altSeqStr);
    altTrans->computePhases();
    if(reverseStrand) altTrans->reverseComplement(altSeqLen);
    Essex::CompositeNode *altTransEssex=altTrans->toEssex();
    altTransEssex->getTag()="mapped-transcript";
    VariantClassifier classifier(variants,VariantClassifier::ALT,*altTrans);
    altTransEssex->append(classifier.makeVariantsNode());
    root->append(altTransEssex);
    ProjectionChecker checker(*refTrans,*altTrans,refSeqStr,refSeq,
			      altSeqStr,altSeq,projectedLab,sensors);
    TranscriptSignals *signals=checker.findBrokenSpliceSites();
    if(!signals) {
      status->prepend("unequal-numbers-of-exons");
      if(!xmlFilename.empty()) writeXML();
      osFBI<<*root<<endl;
      osFBI<<"#===========================================================\n";
      delete altTrans;
      return -1;
    }

    // Enumerate alternative structures
    if(signals->anyBroken()) 
      return enumerateAlts(altTransEssex,signals,status,altSeqStr,altSeqLen,
			   refProtein,refTrans,altTrans,xmlFilename,root,
			   osFBI);

    // Otherwise, projection was successful
    status->prepend("mapped");
    mapped=true;

    // Translate to proteins
    if(refTrans->isCoding()) {
      int oldOrfLen, newOrfLen; float oldStartScore, newStartScore;
      Essex::CompositeNode *upstreamStart=
	orfAnalyzer->earlierStartCodon(*refTrans,refSeqStr,refSeq,
				       *altTrans,altSeqStr,altSeq,
				       *revAlignment,oldOrfLen,newOrfLen,
				       oldStartScore,newStartScore);
      if(upstreamStart) {
	Essex::CompositeNode *changeNode=
	  new Essex::CompositeNode("new-upstream-start-codon");
	changeNode->append("new-start-codon-score",newStartScore);
	changeNode->append("old-start-codon-score",oldStartScore);
	Essex::CompositeNode *lengthNode=
	  new Essex::CompositeNode("ORF-length");
	lengthNode->append(oldOrfLen);
	lengthNode->append("=>");
	lengthNode->append(newOrfLen);
	changeNode->append(lengthNode);
	changeNode->append(upstreamStart);
	status->append(changeNode);
      }

      String refProtein, altProtein;
      checker.translate(*refTrans,*altTrans,refProtein,altProtein);
      
      int ejcDistance;
      switch(nmd.predict(*altTrans,altSeqStr,ejcDistance)) {
      case NMD_NONE: break;
      case NMD_NMD: {
	//status->append("premature-stop","NMD"); break;
	Essex::CompositeNode *stopNode=
	  new Essex::CompositeNode("premature-stop");
	stopNode->append("NMD");
	stopNode->append("EJC-distance",ejcDistance);
	status->append(stopNode);
      } break;
      case NMD_TRUNCATION:  { // ### this is disabled for now
	int matches, len;
	alignProteins(refProtein,altProtein,matches);
	Essex::CompositeNode *fate=new Essex::CompositeNode("premature-stop");
	fate->append("protein-truncation");
	percentMatch(matches,refProtein.length(),altProtein.length(),fate);
	status->append(fate);
      }
	break;
      case NMD_NO_STOP: 
	if(refTrans->hasUTR3()) status->append("nonstop-decay");
	break;
      case NMD_NO_START: status->append("no-start-codon"); break;
      }

      // Check for start codon
      if(startCodonMsg) status->append(startCodonMsg);
      
      // Check for frameshifts
      if(refProtein!=altProtein) {
	int matches, len;
	alignProteins(refProtein,altProtein,matches);
	Essex::CompositeNode *fate=new Essex::CompositeNode("protein-differs");
	percentMatch(matches,refProtein.length(),altProtein.length(),fate);
	status->append(fate);
	Labeling altLab(altSeqLen);
	computeLabeling(*altTrans,altLab);
	checker.checkFrameshifts(projectedLab,altLab,status);
      }
    }
    else { // ref gene is noncoding
      if(!quiet) status->append("noncoding");
      int refOrfLen, altOrfLen;
      Essex::CompositeNode *codingTranscript=
	orfAnalyzer->noncodingToCoding(*refTrans,refSeqStr,refSeq,*altTrans,
				       altSeqStr,altSeq,refOrfLen,altOrfLen);
      if(codingTranscript) {
	Essex::CompositeNode *changeNode=
	  new Essex::CompositeNode("noncoding-to-coding");
	Essex::CompositeNode *lengthNode=
	  new Essex::CompositeNode("ORF-length");
	lengthNode->append(refOrfLen);
	lengthNode->append("=>");
	lengthNode->append(altOrfLen);
	changeNode->append(lengthNode);
	changeNode->append(codingTranscript);
	status->append(changeNode);
      }
    }
    delete altTrans;
  }

  // Flush output>
  if(mapped && status && status->getNumChildren()<2 && quiet) return 0;
  if(!xmlFilename.empty()) writeXML();
  osFBI<<*root<<endl;
  osFBI<<"#===========================================================\n";
  return 0;
}



int FBI::enumerateAlts(Essex::CompositeNode *altTransEssex,
		       TranscriptSignals *signals,
		       Essex::CompositeNode *status,
		       const String &altSeqStr,
		       int altSeqLen,
		       const String &refProtein,
		       const GffTranscript *refTrans,
		       GffTranscript *altTrans,
		       const String &xmlFilename,
		       Essex::CompositeNode *root,
		       ostream &osFBI)
{
  altTransEssex->deleteChild("translation");
  appendBrokenSignals(signals,status);
  EnumerateAltStructures enumerator(*signals,altSeqStr,MAX_SPLICE_SHIFT,
				    MIN_EXON_LEN,MIN_INTRON_LEN,
				    NMD_DISTANCE_PARM,sensors,
				    allowExonSkipping,allowIntronRetention,
				    allowCrypticSites);
  const Vector<AlternativeStructure*> &altStructures=
    enumerator.getAltStructures();
  const int numStruct=altStructures.size();
  if(numStruct>0) {
    status->prepend("splicing-changes");
    Essex::CompositeNode *altStructNode=
      new Essex::CompositeNode("alternate-structures");
    status->append(altStructNode);
    for(Vector<AlternativeStructure*>::const_iterator cur=
	  altStructures.begin(), end=altStructures.end() ; cur!=end ; 
	++cur) {
      const AlternativeStructure &s=**cur;
      Essex::CompositeNode *msg=s.msg;
      s.transcript->loadSequence(altSeqStr);
      s.transcript->computePhases();
      if(reverseStrand) s.transcript->reverseComplement(altSeqLen);
      Essex::CompositeNode *node=s.transcript->toEssex();
      VariantClassifier classifier(variants,VariantClassifier::ALT,
				   *s.transcript);
      node->append(classifier.makeVariantsNode());
      s.reportCrypticSites(node,reverseStrand,altSeqLen);
      if(s.structureChange.anyChange()) {
	Essex::CompositeNode *changeNode=
	  new Essex::CompositeNode("structure-change");
	node->prepend(changeNode);
	if(s.structureChange.crypticSite) 
	  changeNode->append("cryptic-site");
	if(s.structureChange.exonSkipping) 
	  changeNode->append("exon-skipping");
	if(s.structureChange.intronRetention) 
	  changeNode->append("intron-retention");
	if(msg) { changeNode->append(msg); msg=NULL; }
      }
      if(msg) status->append(msg);
      altStructNode->append(node);
      handleProteinFate(s,node,refProtein,altSeqStr,refTrans);
    }
  }
  else status->prepend("no-transcript");
  if(!xmlFilename.empty()) writeXML();
  osFBI<<*root<<endl;
  osFBI<<"#===========================================================\n";
  delete altTrans;
  return 0;
}



void FBI::handleProteinFate(const AlternativeStructure &s,
			    Essex::CompositeNode *node,
			    const String refProtein,
			    const String &altSeqStr,
			    const GffTranscript *refTrans)
{
  switch(s.proteinFate) {
  case NMD_NONE:{      // nothing wrong
    s.transcript->loadSequence(altSeqStr);
    bool identical=s.transcript->getProtein()==refProtein;
    if(identical)
      node->append("fate","identical-protein");
    else {
      int matches, len;
      String altProtein=s.transcript->getProtein();
      alignProteins(refProtein,altProtein,matches);
      Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
      fate->append("protein-differs");
      percentMatch(matches,refProtein.length(),altProtein.length(),
		   fate);
      node->append(fate);
    }
  }break;
  case NMD_NMD: {       // premature stop codon & NMD
    //node->append("fate","NMD");
    Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
    fate->append("NMD");
    fate->append("EJC-distance",s.ejcDistance);
    node->append(fate);
  }break;
  case NMD_TRUNCATION: {// premature stop codon, truncated protein
    int matches, len;
    String altProtein=s.transcript->getProtein();
    alignProteins(refProtein,altProtein,matches);
    Essex::CompositeNode *fate=new Essex::CompositeNode("fate");
    fate->append("protein-truncation");
    percentMatch(matches,refProtein.length(),altProtein.length(),fate);
    node->append(fate);
  }
    break;
  case NMD_NO_STOP:    // no stop codon
    if(refTrans->hasUTR3()) // can't predict if no annotated UTR
      node->append("fate","nonstop-decay");
    break;
  case NMD_NO_START:   // no start codon
    node->append("fate","noncoding");
    break;
  }
}



void FBI::appendBrokenSignals(const TranscriptSignals *signals,
				      Essex::CompositeNode *status)
{
  int numSignals=signals->numSignals();
  for(int i=0 ; i<numSignals ; ++i) {
    TranscriptSignal &signal=(*signals)[i];
    if(!signal.isBroken()) continue;
    String tag;
    SignalType type=signal.getType();
    if(type==GT) tag=signal.weakened ? "weakened-donor" : "broken-donor";
    else if(type==AG) 
      tag=signal.weakened ? "weakened-acceptor" : "broken-acceptor";
    else INTERNAL_ERROR;
    Essex::CompositeNode *node=new Essex::CompositeNode(tag);
    int pos=signal.getPos();
    if(reverseStrand) pos=altSeqLen-pos-1;
    node->append(pos);
    Vector<String> fields; signal.seq.getFields(fields);
    for(Vector<String>::iterator cur=fields.begin(), end=fields.end() ; 
	cur!=end ; ++cur) node->append(*cur);
    //node->append(signal.seq);
    //node->append(signal.refScore);
    //Essex::CompositeNode *cutoffNode=new Essex::CompositeNode("threshold");
    SignalSensor *sensor=sensors.findSensor(type);
    //cutoffNode->append(float(sensor->getCutoff()));
    //node->append(cutoffNode);
    node->append("threshold:"); node->append(float(sensor->getCutoff()));
    status->append(node);
  }
}



void FBI::writeXML()
{
  ofstream os(xmlFilename.c_str());
  root->printXML(os);
  os<<endl;
}



void FBI::append(Essex::CompositeNode *root,const String &tag,
			 const String &message)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(message);
  root->appendChild(node);
}



void FBI::append(Essex::CompositeNode *root,const char *tag,
			 const char *message)
{
  append(root,String(tag),String(message));
}



void FBI::append(Essex::CompositeNode *root,const char *tag,int x)
{
  Essex::CompositeNode *node=new Essex::CompositeNode(tag);
  node->append(x);
  root->appendChild(node);
}



void FBI::append(Essex::CompositeNode *root,const char *tag,
			const String &message)
{
  append(root,tag,message.c_str());
}



void FBI::computeLabeling(GffTranscript &transcript,
				  Labeling &refLab)
{
  const int begin=transcript.getBegin(), end=transcript.getEnd();
  char strand=transcript.getStrand();
  if(strand!='+') throw "only forward-strand features are currently supported";
  int numExons=transcript.getNumExons();
  refLab.asArray().setAllTo(LABEL_INTERGENIC);
  for(int i=begin ; i<end ; ++i) refLab[i]=LABEL_INTRON;
  int phase=0;
  for(int i=0 ; i<numExons ; ++i) {
    GffExon &exon=transcript.getIthExon(i);
    const int begin=exon.getBegin(), end=exon.getEnd();
    for(int j=begin ; j<end ; ++j) {
      refLab[j]=getExonLabel(phase);
      phase=(phase+1)%3;
    }
  }
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    GffExon *UTR=*cur;
    refLab.setIntervalTo(Interval(UTR->getBegin(),UTR->getEnd()),LABEL_UTR);
  }
}



void FBI::mapLabeling(Labeling &from,Labeling &to,
			      const CigarString &cigar)
{
  //  CigarAlignment &align=*cigar.getAlignment();
  CigarAlignment &align=*alignment;
  to.asArray().setAllTo(LABEL_NONE);
  int L=align.length();
  for(int i=0 ; i<L ; ++i) {
    int j=align[i];
    if(j!=CIGAR_UNDEFINED) to[j]=from[i];
  }
  delete &align;
}



bool FBI::mapExon(GffExon &exon,CigarAlignment &align)
{
  int begin=exon.getBegin(), end=exon.getEnd();

  // These two lines map the splice sites across the
  // alignment, then use that to set exon boundaries:
  begin=align.mapApproximate(begin-2,DIR_NONE)+2;
  end=align.mapApproximate(end,DIR_NONE);
  if(begin<0 || end<0) return false;
  exon.setBegin(begin); exon.setEnd(end);
  return true;
}



void FBI::mapTranscript(GffTranscript &refTrans,
				const CigarString &cigar,
				const String &outfile,
				const String &altSeqStr,
				const Sequence &altSeq)
{
  CigarAlignment &align=*alignment;//*cigar.getAlignment();
  Vector<GffExon*> rawExons;
  refTrans.getRawExons(rawExons);
  GffTranscript transcript(refTrans.getTranscriptId(),
			   refTrans.getSubstrate(),
			   refTrans.getStrand(),"FBI");
  transcript.setGeneId(refTrans.getGeneId());
  transcript.setSubstrate(substrate);
  transcript.getSource()="FBI";
  for(Vector<GffExon*>::iterator cur=rawExons.begin(), end=rawExons.end() ;
	cur!=end ; ++cur) {
    GffExon *exon=new GffExon(**cur,transcript);
    if(!mapExon(*exon,align)) INTERNAL_ERROR;
    transcript.addUTR(exon);
  }
  GffTranscript::deleteExons(rawExons);
  if(refTrans.isCoding()) {
    int mappedStartCodon=
      align.mapApproximate(refTrans.getIthExon(0).getBegin(),DIR_LEFT);
    int startCodon=StartCodonFinder::findStartCodon(transcript,
						    transcript.peekUTR(),
						    altSeqStr,
						    mappedStartCodon,
						    sensors);
    if(startCodon!=mappedStartCodon)
      if(startCodon>0) {
	Essex::CompositeNode *newNode
	  =new Essex::CompositeNode("start-codon-change");
	newNode->append("from",mappedStartCodon);
	newNode->append("to",startCodon); 
	startCodonMsg=newNode; }
      else startCodonMsg=new Essex::StringNode("no-start-codon");
    if(startCodon>=0)
      transcript.splitUTRandCDS(altSeqStr,startCodon,sensors.stopCodons);
  }
  transcript.setExonTypes(); transcript.setUTRtypes();
  delete &align;
  ofstream os(outfile.c_str());
  transcript.toGff(os);
}



String FBI::loadSeq(const String &filename)
{
  FastaReader reader(filename);
  String def, seq;
  if(!reader.nextSequence(def,seq)) throw filename+" : cannot read file";
  return seq;
}



String FBI::loadSeq(const String &filename,String &CIGAR)
{
  FastaReader reader(filename);
  String seq, remainder;
  if(!reader.nextSequence(altDefline,seq)) 
    throw filename+" : cannot read file";
  const int L=seq.length();
  FastaReader::parseDefline(altDefline,substrate,remainder);
  if(warningsRegex.search(remainder)) VCFwarnings=warningsRegex[1];
  if(errorsRegex.search(remainder)) VCFerrors=errorsRegex[1];
  if(coordRegex.search(remainder)) globalCoord=coordRegex[1];
  Map<String,String> attr;
  FastaReader::parseAttributes(remainder,attr);
  if(!attr.isDefined("cigar")) 
    throw String("No CIGAR string found on defline: ")+altDefline;
  CIGAR=attr["cigar"];
  parseVariants(attr["variants"],variants,L);
  return seq;
}



GffTranscript *FBI::loadGff(const String &filename)
{
  GffReader reader(filename);
  Vector<GffTranscript*> *transcripts=reader.loadTranscripts();
  const int n=transcripts->size();
  if(n<1) throw filename+" contains no transcripts";
  GffTranscript *transcript=(*transcripts)[0];
  for(int i=1 ; i<n ; ++i) delete (*transcripts)[i];
  delete transcripts;
  transcript->setExonTypes();
  transcript->setUTRtypes();
  return transcript;
}



void FBI::writeProtein(const String &def,const String &protein,
			       const String &filename)
{
  if(filename.empty()) return;
  fastaWriter.writeFasta(def,protein,filename);
}



void FBI::processConfig(const String &filename)
{
  const String path=File::getPath(filename);
  char *oldPath=new char[PATH_MAX];
  getcwd(oldPath,PATH_MAX);
  chdir(path.c_str());

  ConfigFile config(filename);
  MAX_SPLICE_SHIFT=config.getIntOrDie("max-splice-shift");
  MIN_EXON_LEN=config.getIntOrDie("min-exon-length");
  MIN_INTRON_LEN=config.getIntOrDie("min-intron-length");
  NMD_DISTANCE_PARM=config.isDefined("NMD-distance") ?
    config.getIntOrDie("NMD-distance") : 50;
  nmd.setDistParm(NMD_DISTANCE_PARM);

  openPenalty=-config.getFloatOrDie("gap-open-penalty");
  extendPenalty=-config.getFloatOrDie("gap-extend-penalty");
  bandwidth=config.getIntOrDie("bandwidth");
  String matrixFile=config.lookupOrDie("subst-matrix");
  substMatrix=
    new SubstitutionMatrix<float>(matrixFile,AminoAlphabet::global());

  allowExonSkipping=config.getBoolOrDie("allow-exon-skipping");
  allowIntronRetention=config.getBoolOrDie("allow-intron-retention");
  allowCrypticSites=config.getBoolOrDie("allow-cryptic-sites");

  parseConsensusList("donor-consensus",config,sensors.donorConsensuses);
  parseConsensusList("acceptor-consensus",config,sensors.acceptorConsensuses);
  parseConsensusList("stop-codons",config,sensors.stopCodons);
  parseConsensusList("start-codons",config,sensors.startCodons);

  sensors.startCodonSensor=loadModel("start-codon-model",config);
  sensors.shortStartSensor=loadModel("short-start-codon-model",config);
  sensors.stopCodonSensor=loadModel("stop-codon-model",config);
  sensors.donorSensor=loadModel("donor-model",config);
  sensors.acceptorSensor=loadModel("acceptor-model",config);
  sensors.setConsensuses();

  const int MIN_ORF_LEN=config.getIntOrDie("min-orf-length");
  orfAnalyzer=new OrfAnalyzer(sensors,MIN_ORF_LEN);

  chdir(oldPath);
  delete [] oldPath;
}



void FBI::parseConsensusList(const String &tag,ConfigFile &config,
				     Set<String> &into)
{
  String consensusString=config.lookupOrDie(tag);
  Vector<String> fields;
  consensusString.getFields(fields,",");
  for(Vector<String>::const_iterator cur=fields.begin(), end=fields.end() ;
      cur!=end ; ++cur)
    into.insert(*cur);
}



SignalSensor *FBI::loadModel(const String &label,ConfigFile &config)
{
  String filename=config.lookupOrDie(label);
  return SignalSensor::load(filename,garbageCollector);
}



float FBI::alignProteins(const String &refStr,const String &altStr,
			 int &matches)
{
  if(refStr.length()==0 || altStr.length()==0) return 0.0;
  const AminoAlphabet &alphabet=AminoAlphabet::global();
  Sequence refSeq(refStr,alphabet), altSeq(altStr,alphabet);
  BandedSmithWaterman<float> aligner(alphabet,refSeq,altSeq,*substMatrix,
				     openPenalty,extendPenalty,bandwidth);
  Alignment *alignment=aligner.fullAlignment();
  matches=alignment->countMatches();
  float score=float(matches)/refStr.length();
  delete alignment;
  return score;
}



void FBI::percentMatch(int matches,int refLen,int altLen,
		       Essex::CompositeNode *parent)
{
  float percent=float(matches)/max(refLen,altLen);
  Essex::CompositeNode *node=new Essex::CompositeNode("percent-match");
  node->append(String(int(10000*percent+5.0/9)/100.0));
  node->append(String(matches)+"/"+max(refLen,altLen));
  node->append(String("ref=")+refLen);
  node->append(String("alt=")+altLen);
  parent->append(node);
}



void FBI::parseVariants(const String &s,Vector<Variant> &variants,int L)
{
  Vector<String> fields;
  s.getFields(fields,",");
  const int numFields=fields.size();
  for(int i=0 ; i<numFields ; ++i) {
    if(!variantRegex.match(fields[i])) throw "Can't parse variant "+fields[i];
    String id=variantRegex[1], chr=variantRegex[2];
    int refPos=variantRegex[3].asInt(), altPos=variantRegex[4].asInt();
    /*if(reverseStrand) {
      refPos=L-refPos-1;
      altPos=L-altPos-1;
      }*/
    String ref=variantRegex[5], alt=variantRegex[6];
    Variant v(id,chr,refPos,altPos,i);
    v.addAllele(ref); v.addAllele(alt);
    variants.push_back(v);
  }
}



Essex::CompositeNode *FBI::makeEssexVariants()
{
  Essex::CompositeNode *parent=new Essex::CompositeNode("variants");
  const int numVariants=variants.size();
  for(int i=0 ; i<numVariants ; ++i) {
    const Variant &v=variants[i];
    String s=v.id+":"+v.chr+":"+v.refPos+":"+v.altPos+":"+v.alleles[0]
      +":"+v.alleles[1];
    parent->append(s);
  }
  return parent;
}



