#!/usr/bin/perl
#

use strict;
use warnings;

use Data::Dumper;
use RNA;

# 1. Detect MSA format of file 'Rfam.seed'
#
# Note: using FILE_FORMAT_MSA_DEFAULT here checks for all MSA type known to RNAlib
my $msa_type = RNA::file_msa_detect_format("Rfam.seed", RNA::FILE_FORMAT_MSA_DEFAULT);

# 2. Break if file format could not be detected
exit (1) if $msa_type == RNA::FILE_FORMAT_MSA_UNKNOWN;

# 3. open file handle to actually read data from
open my $fh, "<Rfam.seed";

# 4. Create variables to store record data in the following loop
my ($num, $names, $aln, $id, $structure);

# 5. Read records from file as long as file_msa_read_record does not return -1
while ( (($num, $names, $aln, $id, $structure) = RNA::file_msa_read_record($fh, $msa_type)) && ($num != -1)) {

  # skip empty alignments, i.e. number of sequences == 0
  next if $num == 0;

  # print alignment
  print $id, " ", $structure, "\n";

  print Dumper($aln), "\n";
}
