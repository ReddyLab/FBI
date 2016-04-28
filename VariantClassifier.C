/****************************************************************
 VariantClassifier.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantClassifier.H"
using namespace std;
using namespace BOOM;



VariantClassifier::VariantClassifier(const Vector<Variant> &variants,
				     RefAlt refAlt,
				     const GffTranscript &transcript)
{
  classify(variants,refAlt,transcript);
}



const Vector<VariantClassifier::VariantInfo> &VariantClassifier::getCDSvariants() const
{
  return cdsVariants;
}



const Vector<VariantClassifier::VariantInfo> &VariantClassifier::getSpliceSiteVariants() const
{
  return spliceSiteVariants;
}



const Vector<VariantClassifier::VariantInfo> &VariantClassifier::getNearSpliceVariants() const
{
  return nearSpliceVariants;
}



const Vector<VariantClassifier::VariantInfo> &VariantClassifier::getAll() const
{
  return all;
}



void VariantClassifier::classify(const Vector<Variant> &variants,
				 RefAlt refAlt,
				 const GffTranscript &transcript)
{
  int transcriptBegin=transcript.getBegin(), transcriptEnd=transcript.getEnd();
  for(Vector<Variant>::const_iterator cur=variants.begin(), end=
	variants.end() ; cur!=end ; ++cur) {
    const Variant variant=*cur;
    const int pos=refAlt==REF ? variant.refPos : variant.altPos;
    const int endPos=pos+(refAlt==REF ? variant.alleles[0].length() :
			  variant.alleles[1].length());
    Interval variantInterval(pos,endPos);
    VariantInfo info;
    info.variant=variant;
    info.type=getVariantType(variant);
    if(pos<transcriptBegin || pos>=transcriptEnd)
      { info.elem=INTERGENIC; continue; }
    Vector<Interval> introns;
    transcript.getIntrons(introns);
    for(Vector<Interval>::const_iterator cur=introns.begin(), end=
	  introns.end() ; cur!=end ; ++cur) {
      const Interval &intron=*cur;
      int iBegin=intron.getBegin(), iEnd=intron.getEnd();
      Interval ss1(iBegin,iBegin+2), ss2(iEnd-2,iEnd);
      if(variantInterval.overlaps(ss1) || variantInterval.overlaps(ss2))
	{ info.elem=SPLICE_SITE; continue; }
      if(variantInterval.overlaps(intron)) { info.elem=INTRON; continue; }
    }
    for(Vector<GffExon*>::iterator cur=transcript.getExons(), end=
	  transcript.getExonsEnd() ; cur!=end ; ++cur) {
      const GffExon *exon=*cur;
      if(exon->overlaps(variantInterval)) { info.elem=CDS; continue; }
    }
    for(Vector<GffExon*>::iterator cur=transcript.getUTR(), end=
	  transcript.getUTRend() ; cur!=end ; ++cur) {
      const GffExon *exon=*cur;
      if(exon->overlaps(variantInterval)) { info.elem=CDS; continue; }
    }
  }
}


/*
struct VariantInfo {
  Variant variant;
  VariantType type; // SNP/INSERTION/DELETION
  FunctionalElement elem;
  int distanceToSpliceSite;
};
  enum FunctionalElement {
    UTR,
    CDS,
    INTRON,
    SPLICE_SITE,
    INTERGENIC
  };
 */


VariantType VariantClassifier::getVariantType(const Variant &variant) const
{
  const int L1=variant.alleles[0].length(), L2=variant.alleles[1].length();
  if(L1==0) return VARIANT_INSERTION;
  if(L2==0) return VARIANT_DELETION;
  if(L1==1 && L2==1) return VARIANT_SNP;
  return VARIANT_COMPLEX;
}


