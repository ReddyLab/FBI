/****************************************************************
 VariantClassifier.H
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantClassifier_H
#define INCL_VariantClassifier_H
#include <iostream>
#include "Variant.H"

using namespace std;
using namespace BOOM;


struct VariantInfo {
  Variant variant;
  VariantType type; // SNP/INSERTION/DELETION
  FunctionalElement elem;
  int distanceToSpliceSite;
};

class VariantClassifier {
public:
  enum RefAlt {REF, ALT};
  enum FunctionalElement {
    CDS,
    SPLICE_SITE,
    INTRON,
    INTERGENIC
  };
  VariantClassifier(const Vector<Variant> &,RefAlt,const GffTranscript &);
  const Vector<VariantInfo> &getCDSvariants() const;
  const Vector<VariantInfo> &getSpliceSiteVariants() const;
  const Vector<VariantInfo> &getNearSpliceVariants() const;
  const Vector<VariantInfo> &getAll() const;
private:
  Vector<VariantInfo> cdsVaiants, spliceSiteVariants, nearSpliceVariants;
  Vector<VariantInfo> all;
  void classify(const Vector<Variant> &,RefAlt,const GffTranscript &);
};

#endif

