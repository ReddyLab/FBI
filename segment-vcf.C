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
using namespace std;
using namespace BOOM;

struct Variant {
  String chr;
  int begin;
  int end;
  Variant(const String &c,int b,int e) : chr(c), begin(b), end(e) {}
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  Regex gzRegex;
  Vector<Variant> variants;
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("segment-vcf <in.vcf> <binsize> <out.bed>");
  const String &infile=cmd.arg(0);
  const int binSize=cmd.arg(1).asInt();
  const String &outfile=cmd.arg(2);
  const bool isZipped=gzRegex.search(infile);

  // Load all variants from the vcf file
  Vector<String> fields;
  File &f=isZipped ? *new GunzipPipe(infile) : *new File(infile);
  while(!f.eof()) {
    String line=f.getline();
    if(line.empty()) continue;
    if(line[0]=='#') continue;
    line.getFields(fields);
    if(fields.size()<10) continue;
    const String &chr=fields[0];
    int pos=fields[1].asInt()-1;
    const String &ref=fields[3];
    int refLen=ref.length();
    //cout<<chr<<"\t"<<pos<<"\t"<<pos+refLen<<endl;
    variants.push_back(Variant(chr,pos,pos+refLen));
  }

  delete &f;
  return 0;
}

