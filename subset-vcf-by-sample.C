/****************************************************************
 subset-vcf-by-sample.C
 Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  int findIndex(const String &id,const Vector<String> &IDs);
  void emitHeaderLines(const Vector<String> &lines,ostream &);
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
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=0)
    throw String("subset-vcf-by-sample <in.vcf> <sampleID> <out.vcf>");
  const String infile=cmd.arg(0);
  const String wantID=cmd.arg(1);
  const String outfile=cmd.arg(2);
  
  // Open output file
  ofstream os(outfile.c_str());

  // Process input file
  VcfReader reader(infile);
  const Vector<String> &sampleIDs=reader.getSampleIDs();
  int wantIndex=findIndex(wantID,sampleIDs);
  if(wantIndex<0) throw String("Can't find sample ID in VCF file: ")+wantID;
  const Vector<String> &headerLines=reader.getHeaderLines();
  emitHeaderLines(headerLines,os);
  const String &chromLine=reader.getChromLine();
  Variant variant; Vector<Genotype> genotypes;
  while(reader.nextVariant(variant,genotypes)) {
    
  }
  reader.close();

  return 0;
}



void Application::emitHeaderLines(const Vector<String> &lines,ostream &os)
{
  for(Vector<String>::const_iterator cur=lines.begin(), end=lines.end() ;
      cur!=end ; ++cur)
    os<<*cur<<endl;
}



int Application::findIndex(const String &id,const Vector<String> &IDs)
{
  int index=0;
  for(Vector<String>::const_iterator cur=IDs.begin(), end=IDs.end() ;
      cur!=end ; ++cur, ++index) {
    const String &thisOne=*cur;
    if(thisOne==id) return index;
  }
  return -1;
}
