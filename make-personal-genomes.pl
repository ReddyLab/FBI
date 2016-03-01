#!/usr/bin/perl
use strict;
use FileHandle;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use ConfigFile;

die "\n
make-personal-genomes.pl <fbi.config> <vcf-dir> <genes.gff> <out-dir>

Assumptions:
 * VCF files must be zipped with bgzip and have accompanying tbi indexes
 * VCF files must contain chromosome in filename, e.g. chr14, chrX, etc.
 * In VCF files, chrom names don't begin with \"chr\", but in GFF and
     2bit files, they do
 * Your environment variable FBI must point to the FBI dir
 * <out-dir> will be populated with two FASTA files per individual
"
  unless @ARGV==4;
my ($configFile,$vcfDir,$gffFile,$outDir)=@ARGV;

#==============================================================
# First, some initialization
#==============================================================

my $VERBOSE=1;
my $config=new ConfigFile($configFile);
my $CHROM_LENGTHS=$config->lookupOrDie("chr-lengths");
my $TABIX=$config->lookupOrDie("tabix");
my $twoBitFile=$config->lookupOrDie("genome");
my $twoBitToFa=$configFile->lookupOrDie("twoBitToFa");
my $IDfile=$configFile->lookupOrDie("individuals");
my $genderFile=$configFile->lookup("gender");
my %chromLen;
loadChromLengths($CHROM_LENGTHS);
my %chrToVCF;
initChrToVCF($vcfDir);
system("mkdir -p $outDir") unless -e $outDir;
system("rm -f $outDir/errors.txt");
my $MARGIN_AROUND_GENE=1000;
my $refGeneFasta="$outDir/refgene.fasta";
my $altGeneFasta="$outDir/altgene.fasta";
my $tempBedFile="$outDir/temp.bed";
my $geneVcfFile="$outDir/gene.vcf";#"$outDir/gene.vcf.gz";
my $geneGcfFile="$outDir/gene.gcf";#"$outDir/gene.gcf.gz";
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
  my $file1="$outDir/$indiv-1.fasta";
  my $file2="$outDir/$indiv-2.fasta";
  $fastaFiles{$indiv}=[$file1,$file2];
}
my $file1="$outDir/ref-1.fasta";
my $file2="$outDir/ref-2.fasta";
$fastaFiles{"reference"}=[$file1,$file2];

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
  die unless defined($chromLen{$chr});
  if($end>$chromLen{$chr}) { $end=$chromLen{$chr} }
  my $strand=$gene->getStrand();
  my $name=$gene->getId();
  print "$name $chr $begin $end\n";
  writeBed4($chr,$begin,$end,$name,$tempBedFile);
  System("$twoBitToFa -bed=$tempBedFile -noMask $twoBitFile $refGeneFasta");
  my $chrVcfFile=$chrToVCF{$chr};
  writeBed3($chr,$begin,$end,$tempBedFile);
  System("$TABIX -h $chrVcfFile -R $tempBedFile > $geneVcfFile");
  System("$FBI/vcf-to-gcf -i $IDfile -c -v $geneVcfFile $geneGcfFile");
  writeBed6($chr,$begin,$end,$name,$strand,$tempBedFile);
  system("rm $altGeneFasta");
  my $dashY=$genderFile eq "" ? "" : "-y $genderFile";
  System("$FBI/gcf-to-fasta $dashY -r $geneGcfFile $twoBitFile $tempBedFile $altGeneFasta >& $outDir/err.out");
  my $err=`cat $outDir/err.out`;
  if($err=~/error/) { die }
  system("cat $outDir/err.out >> $outDir/errors.txt");
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
      open(FASTA,">>$file") || die $file;
      $def=">$geneID /coord=$coord /margin=$MARGIN_AROUND_GENE /cigar=$cigar";
      $fastaWriter->addToFasta($def,$seq,\*FASTA);
      close(FASTA);
    }
    undef $seq; undef $def;
    undef $indivID; undef $alleleNum ; undef $geneID ; undef $coord;
  }
  $fastaReader->close();

  my $numIndiv=`wc -l $geneGcfFile`-1;
}

#==============================================================
# Clean up
#==============================================================

unlink($refGeneFasta);
unlink($altGeneFasta);
unlink($tempBedFile);
unlink($geneVcfFile);
unlink($geneGcfFile);


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
  if($chr=~/chr(.+)/) { $chr=$1 }
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
  <IN>;
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
#==============================================================
#==============================================================
#==============================================================

