#! /usr/bin/perl
# Converts table of expected snp counters to a distance matrix - p-distance  - between a single sample and all other samples -to allow parallel processing

#This portion is the same for all distances
#*******************************************
use strict;
use Getopt::Long;

my $tableFileRowSample = "";   #file with the table: each row corresponds to a sample, columns correspond to variants

my $row ="";

GetOptions("table=s" => \$tableFileRowSample, "row=s" => \$row);

open ENDWRITER, '>', "${tableFileRowSample}_${row}_DISTANCE_END" or die "error trying to open output file to write an end of a program: $!";

#collect program arguments
if($tableFileRowSample eq "" or $row eq "" )
{
	print "Error: provide the name of expected values table file  --table TABLEFILE and the row of a distance table to be created --row ROW (0:N-1) \n"; 
   	
    exit 1;
}

#open output file for writing a row of distances
if(! open WRITER, '>', "$tableFileRowSample"."ROW_${row}_DISTANCE.csv") 
{
    endProgramError();  	
    die "error trying to open output file: $!";
}


#open input files for reading a table of expected values
my $tableReader; 

if(! open  $tableReader, $tableFileRowSample)
{
    endProgramError();  
    die "error trying to open input table file with expected values: $!";
}

print "\n\n_____________________________\n";
print "Started distance matrix calculation for row $row \n";
print "****************************\n";
#read one line at a time, and compute distances between selected $row and all other samples, to produce $row-th row of the distanc table
my $totalSamples=0;
my $lineNumber=0; #sample number
our @targetRowValues;
our @currentRowValues;
my @distanceMatrixRow; #1d array of distances for this particular row
my $totalPatterns=0;

while (my $line = <$tableReader>)
{
   	
    if($lineNumber == $row)
    {
        chomp $line;
	    @targetRowValues = split(',', $line);
        $totalPatterns =   $#targetRowValues; 	
    }	
    
	$lineNumber++;	
   	
}
close $tableReader;

$totalSamples = $lineNumber;
our $totalAlleles = 2 * $totalPatterns;

if($totalAlleles <=0)
{
    endProgramError();  
    exit (1);
}

if(! open  $tableReader, $tableFileRowSample )
{
    endProgramError();  
    die "error trying to open input table file with expected values: $!";
}

$lineNumber=0;
while (my $line = <$tableReader>)
{
   	
    if($lineNumber != $row)
    {
        chomp $line;
	    @currentRowValues = split(',', $line);
        $distanceMatrixRow[$lineNumber]=getDistance();        
    }
    else
    {
        $distanceMatrixRow[$lineNumber]=0;
    }	
    
	$lineNumber++;		
}
close $tableReader;

#Write distance row into file

for (my $col=0; $col < $totalSamples ; $col ++) 
{
	
        if($col == 0)
        {
            print WRITER $distanceMatrixRow[$col];
        }
        else
        {
            print WRITER ','.$distanceMatrixRow[$col];
        }  
}
close WRITER;


print ENDWRITER 0; 

close ENDWRITER;

print "\n***************************\n";
print "Produced distance row between sample $row and $totalSamples other samples using $totalPatterns variants \n";
print "_____________________________\n\n";

#This one can be changed for each distance type
#************************************ 
sub getDistance
{
	my $ret=0;
    for(my $i=0; $i<=$#currentRowValues; $i++)
    {
        my $first =  int ($targetRowValues[$i]+0.5);
        my $second =int ($currentRowValues[$i]+0.5);
        $ret = $ret +abs ($first - $second);
    }
    
	return $ret/$totalAlleles;
}

sub endProgramError
{
    print ENDWRITER 1; 
}


