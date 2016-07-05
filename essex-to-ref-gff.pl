#!/usr/bin/perl
use strict;
use EssexParser;
use EssexFBI;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my %seen;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  my $fbi=new EssexFBI($root);
  my $transcriptID=$fbi->getTranscriptID();
  next if $seen{$transcriptID};
  $seen{$transcriptID}=1;
  my $status=$fbi->getStatusString();
  #if($status->hasDescendentOrDatum("bad-annotation")) { next }
  #next unless if($status eq "mapped");
  my $transcript=$fbi->getRefTranscript();
  if($transcript) {
    my $id=$transcript->getTranscriptId();
    $transcript->setTranscriptId($id);
    $id=$transcript->getGeneId();
    $transcript->setGeneId($id);
    print OUT $transcript->toGff();
  }
  undef $root; undef $fbi;
}
close(OUT);


