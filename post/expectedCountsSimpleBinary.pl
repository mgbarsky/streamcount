#! /usr/bin/perl
# Converts table of counters into expected values
#writes results into 2 files: one is where rows are samples, and second is where rows are variants
#gets rid of all variants that are not valid
use strict;
use Getopt::Long;

my $tableSamplesInRows = "";   #file with the counts - each row represents a sample, then counts for each k-mer, comma-separated

my $tableSamplesInColumns= ""; #file with the counts -  each row represents a variant and then there counts of this kmer in each sample

GetOptions( "samplesRows=s" => \$tableSamplesInRows,
            "samplesCols=s" => \$tableSamplesInColumns);

#collect program arguments
if($tableSamplesInRows eq "" || $tableSamplesInColumns eq  "" )
{
	print "Error: provide names of all necessary input files: --samplesRows COUNTS CSV FILE BY ROW --sampleCols COUNTS CSV FILE BY COLUMN \n";
    	
    	exit 1;
}


#open output files for writing a new table with expected values
open WRITERBYROW, '>', "$tableSamplesInRows"."EXPECTED.csv" or die "error trying to open output file: $!";
open WRITERBYCOL, '>', "$tableSamplesInColumns"."EXPECTED.csv" or die "error trying to open output file: $!";


my $totalIndividuals;

print "\n".'_____________________________'."\n";
print "Started expected genomic copies calculations by simple binary method\n";
print '_____________________________'."\n";


#***********************************************************************
# 6. create a matrix file where each row is a sample, each column is an expected value of a k-mer's genomic copy
#***********************************************************************
open my $tableReader, $tableSamplesInRows or die "error trying to open input table file: $!";


#for each file and each pattern count: if count=0 put 0, otherwise put 1

my $row=0; #sample id
while (my $line = <$tableReader>)
{
   	chomp $line;
	my @lineArray = split(',', $line);	
	
	for (my $col=0; $col <= $#lineArray; $col ++) 
	{
		my $counter= abs (int ($lineArray[$col]));
        
        my $estCopy = 0;
        if($counter!=0)
        {
            $estCopy=1;
		}

		
        if($col == 0) #first entry on this line
        {
            print WRITERBYROW $estCopy;
        }
        else
        {
		    print WRITERBYROW ",".$estCopy;
        }
	}
    $row++;
	#print "Processed row $row \n";
	print WRITERBYROW "\n";
	print "Written expected for sample $row to the table - ROWS SAMPLES \n";	
}
$totalIndividuals=$row;
close $tableReader;
close WRITERBYROW;

#*******************************
#9. Now we are going to perform the same procedure but now for the transposed table
#*****************************************
#here we can skip invalid rows (removepatterns true - right away)
open my $tableReader, $tableSamplesInColumns or die "error trying to open transposed input table file: $!";
#for each file (col) and each pattern (row) find conditional probabilities CP(gij=0|cij), log CP(gij=1|cij), log CP(gij=2|cij) and general probability P(cij)
#if a valid row (pattern) - compute expected values and add them to the output line-
my $row=0; #kmer id

while (my $line = <$tableReader>)
{
   	
    chomp $line;
    my @lineArray = split(',', $line);    

    my @expectedForThisKmer;

    for (my $col=0; $col <= $#lineArray; $col ++) 
    {
        my $counter= abs (int ($lineArray[$col]));
	    my $estCopy = 0;
        if($counter!=0)
        {
            $estCopy=1;
	    }
        push(@expectedForThisKmer, $estCopy);       
	}

    #determine if variance is non-zero - not all values are the same across all the samples     
       
    
    print WRITERBYCOL $expectedForThisKmer[0];
    for (my $r=1; $r < $totalIndividuals; $r ++) 
    {
        print WRITERBYCOL ",".$expectedForThisKmer[$r];
    }
    print WRITERBYCOL "\n";   
	
    $row++;	
}
close $tableReader;
close WRITERBYCOL;

print "Finished calculating expected values for $row k-mers and $totalIndividuals samples in transposed table\n";

print "\n".'_____________________________'."\n";
print "Finished both tables of expected values with $row total k-mers and $totalIndividuals samples\n";
print '_____________________________'."\n\n";
