/****************************************************************
 segment-vcf.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/Regex.H"
#include "BOOM/VcfReader.H"
#include "BOOM/GffReader.H"
#include "BOOM/Vector.H"
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  bool haveGFF;
  Vector<Interval> features;
  String getChrom(const String &vcfFilename);
  void loadGFF(const String &filename,const String &chr);
  void sortGFF();
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
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
  : gzRegex("\\.gz$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"g:");
  if(cmd.numArgs()!=3)
    throw String("\n\
segment-vcf [options] <in.vcf> <binsize> <out.bed>\n\
    -g file : also respect features given in GFF file\n\
\n\
    * VCF file must be for a single chromosome only\n\
    * VCF file must be sorted by position\n\
");
  const String &infile=cmd.arg(0);
  const int binSize=cmd.arg(1).asInt();
  const String &outfile=cmd.arg(2);
  const bool isZipped=gzRegex.search(infile);

  // Process optional GFF file
  haveGFF=cmd.option('g');
  if(haveGFF) {
    String chr=getChrom(infile);
    loadGFF(cmd.optParm('g'));
    sortGFF();
  }

  // Process VCF file
  VcfReader(infile);
  Variant variant; Vector<Genotype> genotype;
  while(nextVariant(variant,genotype)) {
    
  }

  cout<<"[done]"<<endl;
  return 0;
}



String Application::getChrom(const String &vcfFilename)
{
  VcfReader reader(vcfFilename);
  Variant variant; Vector<Genotype> G;
  if(!reader.nextVariant(variant,G))
    throw "Error loading first variant from VCF file";
  return variant.getChr();
}



void Application::loadGFF(const String &filename,const String &chr)
{
  GffReader reader(filename);
  GffFeature *f;
  while(f=reader.nextFeature()) {
    if(f->getSubstrate()!=chr) continue;
    features.push_back(Interval(f->getBegin(),f->getEnd()));
  }
}



