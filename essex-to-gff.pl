#!/usr/bin/perl
use strict;
use EssexParser;
use EssexFBI;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  my $fbi=new EssexFBI($root);
  my $transcriptID=$fbi->getTranscriptID();
  my $status=$fbi->getStatusString();
  print "$status\n";

#   $transcript=$fbiReport->getMappedTranscript();
#   $statusString=$fbiReport->getStatusString();
#             status=mapped/splicing-changes/no-transcript
#   $bool=$fbiReport->hasBrokenSpliceSite();
#   $array=$fbiReport->getAltTranscripts();

  undef $root; undef $fbi;
}



