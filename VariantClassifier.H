/****************************************************************
 VariantClassifier.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantClassifier_H
#define INCL_VariantClassifier_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/GffTranscript.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;


class VariantClassifier {
public:
  enum RefAlt {REF, ALT};
  enum FunctionalElement {
    UTR,
    CDS,
    INTRON,
    SPLICE_SITE,
    INTERGENIC
  };
  struct VariantInfo {
    Variant variant;
    VariantType type; // SNP/INSERTION/DELETION
    FunctionalElement elem;
    int distanceToSpliceSite;
    VariantInfo() : distanceToSpliceSite(-1) {}
  };
  VariantClassifier(const Vector<Variant> &,RefAlt,const GffTranscript &);
  const Vector<VariantInfo> &getCDSvariants() const;
  const Vector<VariantInfo> &getSpliceSiteVariants() const;
  const Vector<VariantInfo> &getNearSpliceVariants() const;
  const Vector<VariantInfo> &getAll() const;
private:
  Vector<VariantInfo> cdsVariants, spliceSiteVariants, nearSpliceVariants;
  Vector<VariantInfo> all;
  void classify(const Vector<Variant> &,RefAlt,const GffTranscript &);
  VariantType getVariantType(const Variant &) const;
};

#endif

