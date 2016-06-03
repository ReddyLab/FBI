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
  my $fbi=new EssexFBI($root);
  $transcriptID=$fbi->getTranscriptID();
  print "$transcriptID\n";

  undef $root; undef $fbi;
}



