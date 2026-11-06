#!/usr/bin/perl
use strict;
use List::Util qw/max min/;
my $usage =<<USAGE;
example: perl $0 <snp151Common_hg19.txt>
USAGE
die $usage if (@ARGV != 1);

my $snp = shift;
my $last_chr = 'chr';

open (IN, $snp) or die ("can't open the $snp");

open (OUT,">$last_chr.snp") or die ("can't create the file!\n");
while(<IN>){
   chomp;
   my @line = split;
   my ($ALT,$QUAL) = split(',',$line[7],2);
   my $MAF = min split(',',$line[8]);
   if ($line[4] eq '-'){
       print "$ALT\t$QUAL\t\t";
       $ALT =~ tr/ATCGatcg/TAGCtagc/;
       $QUAL =~ tr/ATCGatcg/TAGCtagc/;
       print "$ALT\t$QUAL\n";
   }
   if ($line[0] ne $last_chr){
     close OUT;
     $last_chr = $line[0];
     open (OUT,">$last_chr.snp") or die ("can't create the file!\n");
     print OUT "CHROM\tstart\tend\tID\tREF\tALT\tQUAL\tMAF\n";
   }
   my @p = split(',',$line[8]);
   my $max = max @p;
   my ($sub) = grep {$p[$_] == $max} 0..$#p;
   #print "$sub\t$line[8]\t";
   splice @p,$sub,1;
   my $start = $line[1] + 1;
   my $end = $line[2] + 1;
   print OUT "$line[0]\t$start\t$end\t$line[3]\t$ALT\t$QUAL\t";
   print OUT join ',',@p;
   print OUT "\n";
}
close OUT;

