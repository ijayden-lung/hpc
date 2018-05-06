#!/usr/bin/perl -w

my %retention;
$retention{"exon5"}->{"105,120"} = 100;
$retention{"exon5"}->{"102,190"} = 12;
$retention{"exon5"}->{"178,199"} = 132;
$retention{"exon5"}->{"200,300"} = 12;
$retention{"exon5"}->{"230,300"} = 34;
$retention{"exon5"}->{"240,400"} = 34;



while(my ($exon,$val) = each %retention){
	my %stat;
	my $j = 0;
	foreach my $key (sort keys %$val){
		my ($str,$end) = split /,/,$key;
		if($j == 0){
			$stat{$str} = $end;
		}
		else{
			my $i = 0;
			while(my ($st,$en) = each %stat){
				if($str<=$en && $end>=$st){
					$i++;
					if($end >$en){
						$stat{$st} = $end;
					}
					if($str < $st){
						$stat{$str} = $end;
						delete $stat{$st};
					}
				}
			}
			if($i == 0){
				$stat{$str} = $end;
			}
		}
		$j++;
	}
	while(my ($key,$val) = each %stat){
		print "$key\t$val\n";
	}
}



