#!/usr/bin/perl
use strict;
use FileHandle;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use ConfigFile;

die "\n
make-personal-genomes.pl <fbi.config> <genes.gff> <out-dir>

Assumptions:
 * VCF files must be zipped with bgzip and have accompanying tbi indexes
 * VCF files must contain chromosome in filename, e.g. chr14, chrX, etc.
 * In VCF files, chrom names don't begin with \"chr\", but in GFF and
     2bit files, they do
 * Your environment variable \$FBI must point to the FBI dir
 * <out-dir> will be populated with FASTA files
"
  unless @ARGV==3;
my ($configFile,$gffFile,$outDir)=@ARGV;

#==============================================================
# First, some initialization
#==============================================================


my $DEBUG=0;
my $VERBOSE=1;
my $MARGIN_AROUND_GENE=1000;
my $config=new ConfigFile($configFile);
my $CHROM_LENGTHS=$config->lookupOrDie("chr-lengths");
my $TABIX=$config->lookupOrDie("tabix");
my $twoBitFile=$config->lookupOrDie("genome");
my $twoBitToFa=$config->lookupOrDie("twoBitToFa");
my $IDfile=$config->lookupOrDie("individuals");
my $genderFile=$config->lookup("gender");
my $vcfDir=$config->lookupOrDie("vcf");
my $ploidy=0+$config->lookupOrDie("ploidy");
if($ploidy<1) { die "invalid ploidy" }
my %chromLen;
loadChromLengths($CHROM_LENGTHS);
my %chrToVCF;
initChrToVCF($vcfDir);
system("mkdir -p $outDir") unless -e $outDir;
system("rm -f $outDir/errors.txt");
my $refGeneFasta="$outDir/refgene.fasta";
my $altGeneFasta="$outDir/altgene.fasta";
my $tempBedFile="$outDir/temp.bed";
my $geneVcfFile="$outDir/gene.vcf";#"$outDir/gene.vcf.gz";
my $geneTvfFile="$outDir/gene.tvf";#"$outDir/gene.tvf.gz";
my $FBI=$ENV{"FBI"};
my $fastaWriter=new FastaWriter;

#==============================================================
# Load gene coordinates from GFF file
#==============================================================

my $gffReader=new GffTranscriptReader();
my $genes=$gffReader->loadGenes($gffFile);

#==============================================================
# Make FASTA files for each individual
#==============================================================

my %keepIDs;
loadIDs($IDfile,\%keepIDs);
$keepIDs{"reference"}=1;
my $individuals=getIndividualList($vcfDir);
my $numIndiv=@$individuals;
my %fastaFiles;
for(my $i=0 ; $i<$numIndiv ; ++$i) {
  my $indiv=$individuals->[$i];
  next unless $keepIDs{$indiv};
  $fastaFiles{$indiv}=[];
  for(my $j=1 ; $j<=$ploidy ; ++$j) {
    my $file="$outDir/$indiv-$j.fasta";
    push @{$fastaFiles{$indiv}},$file;
  }
}
$fastaFiles{"reference"}=[];
for(my $j=1 ; $j<=$ploidy ; ++$j) {
  my $file="$outDir/ref-$j.fasta";
  push @{$fastaFiles{"reference"}},$file;
}

#==============================================================
# Process each gene
#==============================================================

my $numGenes=@$genes;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $gene=$genes->[$i];
  my $chr=$gene->getSubstrate();
  my $begin=$gene->getBegin()-$MARGIN_AROUND_GENE;
  my $end=$gene->getEnd()+$MARGIN_AROUND_GENE;
  if($begin<0) { $begin=0 }
  next unless defined($chromLen{$chr});
  if($end>$chromLen{$chr}) { $end=$chromLen{$chr} }
  my $strand=$gene->getStrand();
  my $name=$gene->getId();
  #print "$name $chr $begin $end\n";
  writeBed4($chr,$begin,$end,$name,$tempBedFile);
  System("$twoBitToFa -bed=$tempBedFile -noMask $twoBitFile $refGeneFasta");
  my $chrVcfFile=$chrToVCF{$chr};
  writeBed3($chr,$begin,$end,$tempBedFile);
  System("$TABIX -h $chrVcfFile -R $tempBedFile > $geneVcfFile");
  #System("$FBI/vcf-to-tvf -i $IDfile -c -v $geneVcfFile $geneTvfFile");
  System("$FBI/vcf-to-tvf -i $IDfile -v $geneVcfFile $geneTvfFile");
  writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
  system("rm $altGeneFasta");
  my $dashY=$genderFile eq "" ? "" : "-y $genderFile";
  my $errFile="$outDir/err.out";
  System("$FBI/tvf-to-fasta $dashY -r $geneTvfFile $twoBitFile $tempBedFile $altGeneFasta >& $errFile");
  my $err=`cat $errFile`;
  if($err=~/error/ || $err=~/Abort/) { die $err }
  my (%warnings,%errors);
  loadErrors($errFile,\%warnings,\%errors);
  system("cat $errFile >> $outDir/errors.txt");
  die unless -e $altGeneFasta;
  die if -z $altGeneFasta;
  my $fastaReader=new FastaReader($altGeneFasta);
  while(1) {
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless $def;
    $def=~/>\S+\s+\/individual=(\S+)\s+\/allele=(\d+)\s+\/locus=(\S+)\s+\/coord=(\S+)\s+\/cigar=(\S+)/
      || die "Can't parse defline: $def\n";
    my ($indivID,$alleleNum,$geneID,$coord,$cigar)=($1,$2,$3,$4,$5);
    if($keepIDs{$indivID}) {
      my $file=$fastaFiles{$indivID}->[$alleleNum];
      my $key="$indivID $geneID";
      my $numWarn=0+$warnings{$key};
      my $numErr=0+$errors{$key};
      open(FASTA,">>$file") || die $file;
      $def=">${geneID}_$alleleNum /coord=$coord /margin=$MARGIN_AROUND_GENE /cigar=$cigar /warnings=$numWarn /errors=$numErr";
      $fastaWriter->addToFasta($def,$seq,\*FASTA);
      close(FASTA);
    }
    undef $seq; undef $def;
    undef $indivID; undef $alleleNum ; undef $geneID ; undef $coord;
  }
  $fastaReader->close();
  last if $DEBUG;
}

#==============================================================
# Clean up
#==============================================================

if(!$DEBUG) {
  unlink($refGeneFasta);
  unlink($altGeneFasta);
  unlink($tempBedFile);
  unlink($geneVcfFile);
  unlink($geneTvfFile);
}

#==============================================================
sub getIndividualList {
  my ($vcfDir)=@_;
  my @files=`ls $vcfDir/*.vcf.gz`;
  die "no VCF files found\n" unless @files>0;
  my $file=$files[0];
  chomp $file;
  my $individuals=[];
  open(IN,"cat $file|gunzip|") || die "can't open file $file\n";
  while(<IN>) {
    chomp;
    if(/^\s*#CHROM/) {
      my @fields=split;
      my $numFields=@fields;
      for(my $i=9 ; $i<$numFields ; ++$i)
	{ push @$individuals,$fields[$i] }
      last;
    }
  }
  close(IN);
  return $individuals;
}
#==============================================================
# writeBed4($chr,$begin,$end,$name,$tempBedFile);
sub writeBed4 {
  my ($chr,$begin,$end,$name,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\n";
  close(OUT);
}
#==============================================================
# writeBed3($chr,$begin,$end,$tempBedFile);
sub writeBed3 {
  my ($chr,$begin,$end,$outfile)=@_;
  #if($chr=~/chr(.+)/) { $chr=$1 }
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\n";
  close(OUT);
}
#==============================================================
# writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
sub writeBed6 {
  my ($chr,$begin,$end,$name,$strand,$outfile)=@_;
  open(OUT,">$outfile") || die "Can't write file $outfile\n";
  print OUT "$chr\t$begin\t$end\t$name\t0\t$strand\n";
  close(OUT);
}
#==============================================================
sub System {
  my ($cmd)=@_;
  if($VERBOSE) { print "$cmd\n" }
  system($cmd);
}
#==============================================================
# initChrToVCF($vcfDir);
sub initChrToVCF {
  my ($dir)=@_;
  my @files=`ls $dir/*.vcf.gz`;
  foreach my $file (@files) {
    chomp $file;
    if($file=~/(chr[A-Za-z\d]+)/) {
      my $chr=$1;
      $chrToVCF{$chr}=$file;
    }
  }
}
#==============================================================
sub loadChromLengths
{
  my ($infile)=@_;
  open(IN,$infile) || die $infile;
  #<IN>; no header!
  while(<IN>) {
    chomp;
    my @fields=split;
    next unless @fields>=2;
    my ($chr,$len)=@fields;
    $chromLen{$chr}=$len;
  }
  close(IN);
}
#==============================================================
sub loadIDs
{
  my ($file,$hash)=@_;
  open(IN,$file) || die "can't open file $file\n";
  while(<IN>) {
    chomp;
    next unless $_=~/$\s*(\S+)/;
    my $ID=$1;
    $hash->{$ID}=1;
  }
  close(IN);
}
#==============================================================
sub loadErrors
{
  my ($filename,$warnings,$errors)=@_;
  open(IN,$filename) || die "can't open file: $filename";
  while(<IN>) {
    chomp;
    my @fields=split;
    next unless @fields>=6;
    my ($severity,$type,$indiv,$gene)=@fields;
    my $key="$indiv $gene";
    if($severity=~/WARNING/) { ++$warnings->{$key} }
    elsif($severity=~/ERROR/) { ++$errors->{$key} }
  }
  close(IN);
}
#==============================================================
#==============================================================
#==============================================================

