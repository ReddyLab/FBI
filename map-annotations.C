/****************************************************************
 map-annotation.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/CigarString.H"
#include "BOOM/GffReader.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  int main(int argc,char *argv[]);
private:
  CigarAlignment *alignment;
  void loadAlignment(const String &filename);
  void mapTranscript(const GffTranscript &,ostream &os);
  void mapExon(GffExon &);
};



int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e) 
    { cerr << "STL exception caught in main:\n" << e.what() << endl; }
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}




int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"s:");
  if(cmd.numArgs()!=3)
    throw String("\n\
map-annotations [options] <ref.gff> <ref-to-alt.cigar> <out.gff>\n\
  -s X : change substrate (chrom) to X\n\
\n");
  const String refGffFile=cmd.arg(0);
  const String cigarFile=cmd.arg(1);
  const String outGff=cmd.arg(2);
  bool optS=cmd.option('s');
  String newSubstrate; if(optS) newSubstrate=cmd.optParm('s');

  // Load alignment (CIGAR)
  loadAlignment(cigarFile);

  // Project genes, if any
  ofstream os(outGff.c_str());
  GffReader reader(refGffFile);
  Vector<GffTranscript*> &transcripts=*reader.loadTranscripts();
  for(Vector<GffTranscript*>::iterator cur=transcripts.begin(), end=
	transcripts.end() ; cur!=end ; ++cur) {
    GffTranscript *transcript=*cur;
    if(optS) transcript->setSubstrate(newSubstrate);
    mapTranscript(*transcript,os);
    delete transcript;
  }
  delete &transcripts;

  // Project non-genic elements, if any


  return 0;
}



void Application::loadAlignment(const String &filename)
{
  String line;
  File f(filename);
  line=f.getline();
  f.close();
  CigarString cigar(line);
  alignment=cigar.getAlignment();
}



void Application::mapExon(GffExon &exon)
{
  CigarAlignment &align=*alignment;
  int begin=exon.getBegin(), end=exon.getEnd();

  // These two lines map the splice sites across the
  // alignment, then use that to set exon boundaries:
  begin=align.mapApproximate(begin-2,DIR_NONE)+2;
  end=align.mapApproximate(end+1,DIR_NONE)-1;
  exon.setBegin(begin); exon.setEnd(end);
}



void Application::mapTranscript(const GffTranscript &refTrans,ostream &os)
{
  GffTranscript transcript=refTrans;
  for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	transcript.getExonsEnd() ; cur!=end ; ++cur)
    mapExon(**cur);
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    mapExon(**cur);
  }
  transcript.toGff(os);
}







