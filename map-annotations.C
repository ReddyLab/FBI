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
  GffTranscript *loadGff(const String &filename);
  void mapTranscript(GffTranscript &,const String &cig,const String &outfile);
  void mapExon(GffExon &,CigarAlignment &);
};



int main(int argc,char *argv[])
  {
    try
      {
	Application app;
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




int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"s:");
  if(cmd.numArgs()!=3)
    throw String("\n\
map-annotations [options] <ref.gff> <ref-to-alt.cigar> <out.gff>\n\
\n\
  -s X : change substrate (chrom) to X\n\
\n\
  NOTE: ref.gff must contain only one transcript\n\
\n");
  const String refGffFile=cmd.arg(0);
  const String cigarFile=cmd.arg(1);
  const String outGff=cmd.arg(2);

  // Load reference GFF and CIGAR string
  GffTranscript *refTrans=loadGff(refGffFile);
  if(cmd.option('s')) refTrans->setSubstrate(cmd.optParm('s'));
  String CIGAR;
  File cigar(cigarFile);
  CIGAR=cigar.getline();
  cigar.close();

  // Project the reference GFF over to an alternate GFF
  mapTranscript(*refTrans,CIGAR,outGff);

  return 0;
}



void Application::mapExon(GffExon &exon,CigarAlignment &align)
{
  int begin=exon.getBegin(), end=exon.getEnd();

  // These two lines map the splice sites across the
  // alignment, then use that to set exon boundaries:
  begin=align.mapApproximate(begin-2,DIR_NONE)+2;
  end=align.mapApproximate(end+1,DIR_NONE)-1;
  exon.setBegin(begin); exon.setEnd(end);
}



void Application::mapTranscript(GffTranscript &refTrans,const String &cig,
				 const String &outfile)
{
  GffTranscript transcript=refTrans;
  CigarString cigar(cig);
  CigarAlignment &align=*cigar.getAlignment();
  for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	transcript.getExonsEnd() ; cur!=end ; ++cur)
    mapExon(**cur,align);
  for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	transcript.getUTRend() ; cur!=end ; ++cur) {
    mapExon(**cur,align);
  }
  delete &align;
  ofstream os(outfile.c_str());
  transcript.toGff(os);
}



GffTranscript *Application::loadGff(const String &filename)
{
  GffReader reader(filename);
  Vector<GffTranscript*> *transcripts=reader.loadTranscripts();
  const int n=transcripts->size();
  if(n<1) throw filename+" contains no transcripts";
  GffTranscript *transcript=(*transcripts)[0];
  for(int i=1 ; i<n ; ++i) delete (*transcripts)[i];
  delete transcripts;
  return transcript;
}



