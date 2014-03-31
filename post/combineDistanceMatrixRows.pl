#! /usr/bin/perl
# Converts rows of distances into a combined table
#to use in phylogenetics programs - modify according to : http://evolution.genetics.washington.edu/phylip/doc/distance.html

#This portion is the same for all distances
#*******************************************
use strict;
use Getopt::Long;

my $tableFileRowSample = "";   #file with the table: each row corresponds to a sample, columns correspond to variants

GetOptions("table=s" => \$tableFileRowSample);

#collect program arguments
if($tableFileRowSample eq "" )
{
	print "Error: provide the name of expected values table file  --table TABLEFILE \n";    	
    exit 1;
}

#open output file for writing a new table of distances in 
my $outputFileName = "$tableFileRowSample"."DISTANCE_NO_HEADER_P-distance.csv";
open WRITER, '>', $outputFileName or die "error trying to open output file $outputFileName: $!";

#open first input row file to know how many rows of distance matrix to consider (as many rows as columns) 

my $rowReader; 

open my $rowReader, "$tableFileRowSample"."ROW_0_DISTANCE.csv" or die "error trying to open input file with distance for row 0: $!";
print "\n\n_____________________________\n";
print "Started combining rows of distance matrix\n";
print "****************************\n";
#read one line at a time, and compute distances between selected $row and all other samples, to produce $row-th row of the distanc table
my $totalSamples=0;

while (my $line = <$rowReader>)
{    
    chomp $line;
    my @lineArray = split(',', $line);
    $totalSamples =   $#lineArray+1;	
}
close $rowReader;

#now in a loop from 0 to $totalSamples-1 add precomputed distance rows to the output table

for(my $row=0; $row < $totalSamples; $row++)
{
    my $tableName = "$tableFileRowSample"."ROW_${row}_DISTANCE.csv";
    open my $rowReader, "$tableFileRowSample"."ROW_${row}_DISTANCE.csv" or die "error trying to open input file $tableName with distance for row $row : $!";
    my $line = <$rowReader>;
    print WRITER $line."\n";
    close $rowReader;
}


print "\n***************************\n";
print "Produced distance table for $totalSamples samples by combining separate distance matrix rows \n";
print "_____________________________\n\n";


