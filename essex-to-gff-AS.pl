#!/usr/bin/perl
use strict;
use EssexParser;
use EssexFBI;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff> <allele#>\n" unless @ARGV==3;
my ($infile,$outfile,$hap)=@ARGV;

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
  #if($status eq "mapped") {
  my $transcript=$fbi->getMappedTranscript();
  if($transcript) {
    my $id=$transcript->getTranscriptId();
    $id.="_$hap";
    $transcript->setTranscriptId($id);
    $id=$transcript->getGeneId();
    $id.="_$hap";
    $transcript->setGeneId($id);
    print OUT $transcript->toGff();
  }
  if($status eq "splicing-changes") {
    my $transcripts=$fbi->getAltTranscripts();
    my $n=@$transcripts;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $transcript=$transcripts->[$i];
      my $id=$transcript->getTranscriptId();
      if($transcript->getSource() eq "SIMULATION") { $id="SIM$i\_$id\_$hap" }
      else { $id="ALT$i\_$id\_$hap" }
      $transcript->setTranscriptId($id);
      $id=$transcript->getGeneId();
      $id.="_$hap";
      $transcript->setGeneId($id);
      print OUT $transcript->toGff();
    }
  }
  elsif($status eq "mapped") {
    my $newStartTranscript=
      $root->pathQuery("report/status/new-upstream-start-codon/transcript");
    if($newStartTranscript) {
      my $id=$transcript->getTranscriptId();
      $id="NEWSTART_$id\_$hap";
      $transcript->setTranscriptId($id);
      $id=$transcript->getGeneId();
      $id.="_$hap";
      $transcript->setGeneId($id);
      print OUT $transcript->toGff();
    }
  }
  undef $root; undef $fbi;
}
close(OUT);
print "[done]\n";

