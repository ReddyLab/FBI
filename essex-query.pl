#!/usr/bin/perl
use strict;
use EssexParser;
use ProgramName;

# e.g., "report/status"
#    operators: >,>=,<,<=,=,!=,~
#    "~" means "contains substring"
#    i.e., probe/sequence~CCTAGCAGT

my $name=ProgramName::get();
die "$name <in.essex> <query>\n" unless @ARGV==2;
my ($infile,$query)=@ARGV;

my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  if(EssexNode::isaNode($root)) {
    my $array=$root->findDescendents($query);
    foreach my $elem (@$array) {
      $elem->print(\*STDOUT);
      print "\n";
    }
  }
}
