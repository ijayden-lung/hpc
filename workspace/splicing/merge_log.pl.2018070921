#!/usr/bin/perl -w
#
my ($dir,$exp) = @ARGV;
if(!defined $dir){
	$dir = '.';
}
if(!defined $exp){
	$exp = 'STAR';
}
my %hash;
my @sample;
my %hash2;
foreach my $file (glob "$dir/STAR/*/Log.final.out"){
	my $sam = (split /\//,$file)[-2];
	push @sample,$sam;
	open FILE,$file;
	my $i=0;
	while(<FILE>){
		chomp;
		$i++;
		my ($row,$num)  = split /\|/;
		if(defined $num){
			$row =~ s/\s+//;
			$num =~ s/\s+//;
		}
		else{
			$row = '';
			$num = '';
		}
		
		$hash{$row}->{$sam} = $num;
		$hash2{$row} = $i;
	}
}

open OUT,">$dir/Stat_$exp.tsv";
print OUT "\t";
print OUT join "\t", @sample;
print OUT "\n";
foreach my $row (sort{$hash2{$a}<=>$hash2{$b}} keys %hash2){
	my $val = $hash{$row};
	print OUT "$row";
	foreach my $sam (@sample){
		print OUT "\t$val->{$sam}";
	}
	print OUT "\n";
}
