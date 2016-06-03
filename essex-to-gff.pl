#!/usr/bin/perl
use strict;
use EssexParser;
use EssexFBI;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <in.essex> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  my $fbi=new EssexFBI($root);
  my $transcriptID=$fbi->getTranscriptID();
  my $status=$fbi->getStatusString();
  if($status eq "mapped") {
    my $transcript=$fbi->getMappedTranscript();
    print OUT $transcript->toGff();
  }
  elsif($status eq "splicing-changes") {
    my $transcripts=$fbi->getAltTranscripts();
    my $n=@$transcripts;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $transcript=$transcripts->[$i];
      my $id=$transcript->getTranscriptId();
      $id="ALT$i\_$id";
      $transcript->setTranscriptId($id);
      print OUT $transcript->toGff();
    }
  }
  undef $root; undef $fbi;
}
close(OUT);


