/****************************************************************
 Variant.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Variant_H
#define INCL_Variant_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Comparator.H"
#include "BOOM/Vector.H"
using namespace std;
using namespace BOOM;



struct Variant {
  String id, chr;
  Vector<String> alleles; // alleles[0] is the ref; multiple alts are OK
  int pos, i;
  Variant(const String &id,const String &chr,int pos,int i=0)
    : id(id), chr(chr), pos(pos), i(i) {}
  void addAllele(const String &a) { alleles.push_back(a); }
  void printOn(ostream &os) const {
    os<<id<<":"<<chr<<":"<<pos;
    for(Vector<String>::const_iterator cur=alleles.begin(), end=alleles.end() ;
	cur!=end ; ++cur) { String s=*cur; if(s=="") s=".";os<<":"<<s; }
  }
  void printAllele(int allele,ostream &os) const {
    os<<id<<":"<<chr<<":"<<pos<<":"<<(alleles[0]==""?".":alleles[0].c_str())
      <<":"<<(alleles[allele]==""?".":alleles[allele].c_str()); }
  bool identical(const Variant &other) {
    if(chr!=other.chr || pos!=other.pos || alleles[0]!=other.alleles[0] ||
       alleles.size()!=other.alleles.size())
      return false;
    DirectComparator<String> cmp;
    Vector<String> A=alleles, B=other.alleles; A.shift(); B.shift();
    VectorSorter<String> As(A,cmp), Bs(B,cmp); 
    As.sortAscendInPlace(); Bs.sortAscendInPlace();
    const int N=A.size();
    for(int i=0 ; i<N ; ++i) if(A[i]!=B[i]) return false;
    return true;
  }
  bool identical(int state,const Variant &other,int otherState) {
    return chr==other.chr && pos==other.pos && alleles[0]==other.alleles[0]
      && alleles[state]==other.alleles[otherState]; }
  int refLen() const { return alleles[0].length(); }
  int refEnd() const { return pos+alleles[0].length(); }
  bool overlaps(const Variant &other) const { 
    const int end=pos+alleles[0].length();
    const int otherEnd=other.pos+other.alleles[0].length();
    return pos<otherEnd && other.pos<end;
  }
  bool covers(const Variant &other) const {
    const int end=pos+alleles[0].length();
    const int otherEnd=other.pos+other.alleles[0].length();
    return pos<=other.pos && end>=otherEnd;
  }
  /*bool deletion(int allele) const {
    return alleles[allele].length()<alleles[0].length(); }*/
  bool insertion(int allele) const {
    return alleles[0].length()==0; } // "pure" insertion
    //return alleles[allele].length()>alleles[0].length(); }
  bool SNP(int allele) const {
    return alleles[0].length()==1 && alleles[allele].length()==1; }
  void trim();
};
ostream &operator<<(ostream &os,const Variant &v) { v.printOn(os); return os; }


struct VariantComp : Comparator<Variant> { // For sorting variants by position
  bool equal(Variant &a,Variant &b)
  { return a.chr==b.chr && a.pos==b.pos && a.refLen()==b.refLen(); }
  bool greater(Variant &a,Variant &b) 
  { return a.chr>b.chr || a.chr==b.chr && 
      (a.pos>b.pos || a.pos==b.pos && a.refLen()>b.refLen()); }
  bool less(Variant &a,Variant &b)    
  { return a.chr<b.chr || a.chr==b.chr && 
      (a.pos<b.pos || a.pos==b.pos && a.refLen()<b.refLen()); }
};




#endif

