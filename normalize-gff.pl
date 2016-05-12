#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name infile.gff > outfile.gff\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($infile);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i) {
  my $transcript=$transcripts->[$i];
  my $gff=$transcript->toGff();
  my $array=$transcript->parseExtraFields();
  my $hash=$transcript->hashExtraFields($array);
  my $type=$hash->{"type"};
  if(!$type) { $type=$hash->{"gene_type"} }
  if($type && $type ne "protein_coding") { $transcript->becomeNoncoding() }
#   $transcript->setExtraFieldsFromKeyValuePairs(\@array); # [key,value]
#   $transcript->setExtraFields($string);
  print "$gff";
}



