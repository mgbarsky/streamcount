#! /usr/bin/perl
# Gives a table of counts per each pattern line
use strict;
use Getopt::Long;

my $countsDir = "";
my $fileListFile = "";   #text file with list of all files to process
my $fileExtension = "";   #table of mappings from each lineid, posid -> countID, rcCountID - positions in the array of counts, (repeated) 


GetOptions("folder=s" => \$countsDir,
            "files=s" => \$fileListFile,
           "ext=s" => \$fileExtension,
            );

#collect program arguments
if($fileListFile eq "" )
{
	print "Error: provide the name of file with file names\n";    	
    	exit 1;
}

#open output file for writing
open WRITER, '>', "$countsDir".'/'."TABLE_LAMBDAS.csv" or die "error trying to open output file: $!";



#read each counters file, find mode after min and write it to the output file, one mode per line

open my $linesReader, $fileListFile or die "error trying to open file list: $!";
my $row = 0;

while (my $line = <$linesReader>)
{
   	chomp $line;

    my @counters;
    my @distribution;

    my $maxCount=0;
	#we now know the file name where the counters are stored
	my $fileName = $countsDir.'/'.$line.$fileExtension;
	open my $fileReader, $fileName or die "error trying to open file $fileName with counts: $!";	
	
	while(my $countline = <$fileReader>)
	{
		chomp $countline;
        my $current_count = int $countline;
        push (@counters,$current_count); 
        if(! $distribution[$current_count])
        {
            $distribution[$current_count]=1;
        }
        else
        {
            $distribution[$current_count]++;
        }
        if($current_count> $maxCount)
        {
            $maxCount = $current_count;
        }
	}

	close $fileReader;
	
	
    #now go through distribution and find the first min - to cut off all the pick (error counts) around zero and one
    my $minFound = 0;
    my $cutoffCount=0;
    for(my $c=1; $c<=$maxCount-1 and !$minFound; $c++)
    {
        my $prevBucketTotal =0; 
        if($distribution[$c])
        {
            $prevBucketTotal=$distribution[$c];
        }
        my  $currBucketTotal =0;
        if( $distribution[$c+1])
        {
             $currBucketTotal=$distribution[$c+1];
        }
        if($currBucketTotal > $prevBucketTotal)
        {
            $cutoffCount=$c+1;
            $minFound=1;
        }
    }
    
    my $maxBucketSeenSoFar=0;
    my $mode=0;
    #find the largest bucket after cutoff
    for(my $c=$cutoffCount; $c<=$maxCount; $c++)
    {
        my $currBucketTotal = $distribution[$c];
        if($currBucketTotal > $maxBucketSeenSoFar)
        {
            $maxBucketSeenSoFar=$currBucketTotal;
            $mode=$c;
        }
    }
    
	print WRITER $mode."\n";
	
	$row++;
	
}
close $linesReader;

print STDERR "Processed $row samples, mode is written to file ".  "$countsDir".'/TABLE_LAMBDAS.csv'."\n";


