#!/usr/bin/perl 

use strict;
use warnings;

unless($#ARGV==1){
	print "USAGE: perl parse_pdb.pl file.pdb chainID\n";
	exit;
}

my $id = $ARGV[0];
my $chain = $ARGV[1];
open(IN, $id);
my @pdb;
my @seq;
while(<IN>){
	my @col = split;
	next unless $col[0] eq 'ATOM' and $col[2] eq 'CA' and $col[4] eq $chain;
	push @pdb, $col[3];
}
&sequence(@pdb);
&wt_eco(@pdb);

#retrieving fasta chain from pdb file
sub sequence{
my @pdb1 = @_;
foreach my $res(@pdb1){
	$res =~ s/ALA/A/;
	$res =~ s/ARG/R/;
	$res =~ s/ASN/N/;
	$res =~ s/ASP/D/;
	$res =~ s/CYS/C/;
	$res =~ s/GLU/E/;
	$res =~ s/GLN/Q/;
	$res =~ s/GLY/G/;
	$res =~ s/HIS/H/;
	$res =~ s/ILE/I/;
	$res =~ s/LEU/L/;
	$res =~ s/LYS/K/;
	$res =~ s/MET/M/;
	$res =~ s/PHE/F/;
	$res =~ s/PRO/P/;
	$res =~ s/SER/S/;
	$res =~ s/THR/T/;
	$res =~ s/TRP/W/;
	$res =~ s/TYR/Y/;
	$res =~ s/VAL/V/;
	push(@seq, $res);
}
my $seq = join('', @seq);
my $fasta = ">".$id."\n".$seq."\n";
print "\n\n", $fasta,"\n\n";
my $len = length($seq);
print "Number of Amino acids: $len\n\n";
}

#calculate molecular weight and extinction coeeficient
sub wt_eco{
	my @pdb2 = @_;
	my $tyr = 0;
	my $trp = 0;
	my $cys = 0;

	foreach my $res1(@pdb2){
	if ($res1 eq 'TYR'){$tyr++;}
	if ($res1 eq 'TRP'){$trp++;}
	if ($res1 eq 'CYS'){$cys++;}
	}
	my $e_pro = ($tyr * 1490)+($trp * 5500)+($cys * 125);

	my @wt;
	foreach my $res2(@pdb2){

	$res2 =~ s/ALA/71.078/;
        $res2 =~ s/ARG/156.187/;
        $res2 =~ s/ASN/114.103/;
        $res2 =~ s/ASP/115.088/;
        $res2 =~ s/CYS/103.144/;
        $res2 =~ s/GLU/129.115/;
        $res2 =~ s/GLN/128.130/;
        $res2 =~ s/GLY/57.052/;
        $res2 =~ s/HIS/137.141/;
        $res2 =~ s/ILE/113.159/;
        $res2 =~ s/LEU/113.159/;
        $res2 =~ s/LYS/128.174/;
        $res2 =~ s/MET/131.198/;
        $res2 =~ s/PHE/147.176/;
        $res2 =~ s/PRO/97.116/;
        $res2 =~ s/SER/87.078/;
        $res2 =~ s/THR/101.105/;
        $res2 =~ s/TRP/186.213/;
        $res2 =~ s/TYR/163.176/;
        $res2 =~ s/VAL/99.132/;
	push(@wt, $res2);

}
chomp @wt;
my $sum = 0;
$sum += $_ for(@wt);
print "Molecular weight: $sum", "\n\n";
my $eco = sprintf "%0.3f", ($e_pro/$sum);
print "Extinction Coefficient: $eco/m.cm\n\n";
}

