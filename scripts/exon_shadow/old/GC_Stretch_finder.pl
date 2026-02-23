#!/usr/bin/perl
use strict;
use warnings;

# Check for proper command-line usage
if (@ARGV != 1) {
    die "Usage: $0 <input_file.fasta>\n";
}

my $input_file = $ARGV[0];

# Open the multifasta file for reading
open my $fasta_fh, '<', $input_file or die "Cannot open file $input_file: $!";

# Initialize variables
my $current_header = "";
my $current_sequence = "";
my %species_lengths;  # Hash to store lengths for each species

# Process the multifasta file
while (my $line = <$fasta_fh>) {
    chomp $line;

    if ($line =~ /^>/) {
        # If it's a header line, process the previous sequence
        process_sequence($current_header, $current_sequence, \%species_lengths) if $current_header;

        # Start a new sequence
        $current_header = $line;
        $current_sequence = "";
    } else {
        # Concatenate sequence lines
        $current_sequence .= $line;
    }
}

# Process the last sequence in the file
process_sequence($current_header, $current_sequence, \%species_lengths) if $current_header;

# Close the file handle
close $fasta_fh;

# Print the average length for each species
foreach my $species (keys %species_lengths) {
    my $average_length = $species_lengths{$species}->{'total_length'} / $species_lengths{$species}->{'sequence_count'};
    print "$species $average_length\n";
}

sub process_sequence {
    my ($header, $sequence, $species_lengths) = @_;

    # Extract species name from the header (assuming header format ">Species_name ...")
    my ($species) = $header =~ /^>(\S+)/;

    # Find and calculate G/C stretches
    my @gc_stretches = $sequence =~ /([GC]{3,})/ig;
    foreach my $gc_stretch (@gc_stretches) {
        $species_lengths->{$species}->{'total_length'} += length($gc_stretch);
        $species_lengths->{$species}->{'sequence_count'}++;
    }
}
