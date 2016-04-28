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

enum FunctionalElement {
  FE_CDS,
  FE_SPLICE_SITE,
  FE_NEAR_SPLICE_SITE,
  FE_INTRON,
  FE_INTERGENIC
};

struct VariantInfo {
  Variant variant;
  VariantType type; // SNP/INSERTION/DELETION
  FunctionalElement elem;
};

class VariantClassifier {
public:
  enum RefAlt {REF, ALT};
  VariantClassifier(const Vector<Variant> &,RefAlt,const GffTranscript &);
  const Vector<VariantInfo> &getCDSvariants() const;
  const Vector<VariantInfo> &getSpliceSiteVariants() const;
  const Vector<VariantInfo> &getNearSpliceVariants() const;
  const Vector<VariantInfo> &getAll() const;
private:
  Vector<VariantInfo> cdsVaiants, spliceSiteVariants, nearSpliceVariants;
  Vector<VariantInfo> all;
};

#endif

