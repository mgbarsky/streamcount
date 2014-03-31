#! /usr/bin/perl
# Gives a table of counts per each sample - input - counts files - already with rc
use strict;
use Getopt::Long;

my $fileListFile = "";   #file with the names of the files
my $fileExtension = "";   #extension suffix to be added to each file name
my $inputFolder=""; #input folder to be added as a prefix to each file name
my $outputFileName="";

GetOptions("files=s" => \$fileListFile,
            "ext=s" => \$fileExtension,
           "folder=s" => \$inputFolder,
            "output=s" => \$outputFileName);

#collect program arguments
if($fileListFile eq "" or $outputFileName eq "" )
{    
	print "Error: provide names of file list file and the output file name\n";
    	
    	exit 1;
}


#open output file for writing a table
open TABLEWRITER, '>', $outputFileName or die "error trying to open output file: $!";

my $countsReader; #file handle to read each counts file

#read all file names into an array of files
my @files; #holds all file names provided in  @fileList file
open my $filesReader, $fileListFile or die "error trying to open file list file: $!";
while (my $fileName = <$filesReader>)
{
   	chomp $fileName;
	push (@files, $fileName);   # we will add it to the array.		
}

close $filesReader;



my $i; #position in an array of patterns
#now we are going in the loop and for each file generate a line with all counters
for (my $fn=0; $fn<= $#files; $fn++)
{	
	#we read counters from each counters file, and as we read - we write it as a comma-separated line in a new table file
    my $COUNTERS_FILE= "${inputFolder}/".$files[$fn].$fileExtension;
	open my $countsReader, $COUNTERS_FILE or die "error trying to open counters file $COUNTERS_FILE: $!";
	my $lineID=0;
	while (my $counter = <$countsReader>)
	{
	   	chomp $counter;
		if($lineID == 0)
        {
            print TABLEWRITER $counter;
        }
        else
        {
            print TABLEWRITER ','.$counter;
        }
		$lineID++;	
	}
	print TABLEWRITER "\n";
	close $countsReader;

	
	print "Processed file $fn with $lineID counts \n";	
}



