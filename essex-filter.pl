#!/usr/bin/perl
use strict;
use EssexParser;
use ProgramName;
$|=1;

my $slurm=$ENV{"SLURM_JOB_ID"};
print "SLURM $slurm\n";
system("hostname");

my $name=ProgramName::get();
die "$name <in.essex> <out.essex>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write file: $outfile\n";
my $badAnnos=0; my $vcfErrors=0; my $mapped=0; my $kept=0;
my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $status=$report->findChild("status");
  next unless $status;
  if($status->hasDescendentOrDatum("bad-annotation"))
     { ++$badAnnos; next }
  if($status->hasDescendentOrDatum("too-many-vcf-errors"))
     { ++$vcfErrors; next }
  if($status->hasDescendentOrDatum("mapped"))
     { ++$mapped; } #next }
  $report->print(\*OUT);
  print OUT "\n";
  ++$kept;
}
close(OUT);

print "bad-annotation $badAnnos\n";
print "too-many-vcf-errors $vcfErrors\n";
print "mapped $mapped\n";
print "retained $kept\n";

