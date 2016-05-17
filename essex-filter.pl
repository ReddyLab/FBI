#!/usr/bin/perl
use strict;
use EssexParser;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.essex> <out.essex>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write file: $outfile\n";
my $badAnnos=0; my $vcfErrors=0; my $mapped=0;
my $parser=new EssexParser($infile);
while(1) {
  my $report=$parser->nextElem();
  last unless $report;
  my $status=$report->findChild("status");
  next unless $status;
  my $code=$status->getIthElem(0);
  next unless $code;
  if($code eq "bad-annotation") { ++$badAnnos; next }
  if($code eq "too-many-vcf-errors") { ++$vcfErrors; next }
  if($code eq "mapped") { ++$mapped; next }
  $report->print(\*OUT);
  print OUT "\n";
}
close(OUT);

print "bad-annotation $badAnnos\n";
print "too-many-vcf-errors $vcfErrors\n";
print "mapped $mapped\n";

