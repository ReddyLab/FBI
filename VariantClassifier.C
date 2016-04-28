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



const Vector<VariantInfo> &VariantClassifier::getCDSvariants() const
{
  return cdsVariants;
}



const Vector<VariantInfo> &VariantClassifier::getSpliceSiteVariants() const
{
  return spliceSiteVariants;
}



const Vector<VariantInfo> &VariantClassifier::getNearSpliceVariants() const
{
  return nearSpliceVariants;
}



const Vector<VariantInfo> &VariantClassifier::getAll() const
{
  return alt;
}



void VariantClassifier::classify(const Vector<Variant> &variants,
				 RefAlt refAlt,
				 const GffTranscript &transcript)
{
  
}




