#! /usr/bin/perl
# Converts table of counters into expected values
#writes results into 2 files: one is where rows are samples, and second is where rows are variants
#gets rid of all variants that are not valid
use strict;
use Getopt::Long;

my $tableSamplesInRows = "";   #file with the counts - each row represents a sample, then counts for each k-mer, comma-separated
my $lambdaFile = ""; #file with lambda values (in the first column) per file
my $tableSamplesInColumns= ""; #file with the counts -  each row represents a variant and then there counts of this kmer in each sample

GetOptions( "samplesRows=s" => \$tableSamplesInRows,
            "samplesCols=s" => \$tableSamplesInColumns,
	        "lambda=s" => \$lambdaFile);

#collect program arguments
if($tableSamplesInRows eq "" || $tableSamplesInColumns eq  "" || $lambdaFile eq  "")
{
	print "Error: provide names of all necessary input files: --samplesRows COUNTS CSV FILE BY ROW --sampleCols COUNTS CSV FILE BY COLUMN --lambda FILE WITH MEAN PER SAMPLE\n";
    	
    	exit 1;
}


#open output files for writing a new table with expected values
open WRITERBYROW, '>', "$tableSamplesInRows"."EXPECTED.csv" or die "error trying to open output file: $!";
open WRITERBYCOL, '>', "$tableSamplesInColumns"."EXPECTED.csv" or die "error trying to open output file: $!";

open TEMPWRITER, '>', "$tableSamplesInRows"."TEMP" or die "error trying to open output file: $!";

open REMOVEWRITER, '>', "$tableSamplesInRows"."REMOVE" or die "error trying to open output file: $!";

#***********************
#1. Array of lambdas for each sample
#***************************
my @lambdaPerFile;
my $totalIndividuals = 0;
#read lambda values into an array: this array is only $totalIndividuals long

my $min = 100000;
my $max = 0;
open my $lambdaPerFileReader, $lambdaFile or die "error trying to open input file with lambda per file: $!";
while (my $line = <$lambdaPerFileReader>)
{
   	chomp $line;
	my @lineArray = split(',', $line);
	
    my $mean = $lineArray[0];

    if($mean > $max)
    {
        $max=$mean;
    }
    if($mean < $min)
    {
        $min=$mean;
    }
	push (@lambdaPerFile, $mean);	#pushing average
	$totalIndividuals++;		
}
print  "There are total $totalIndividuals samples. Mode per sample ranges from $min to $max\n";
close $lambdaPerFileReader;

#****************************
# 2. Total count of each k-mer across all samples. We are going to use the tableSamplesInColumns
# In this table - each k-mer is in a separate line, so we just add up the values in each line
#******************************
#open input files for reading a table of counters
my $tableReader; 
open $tableReader, $tableSamplesInColumns or die "error trying to open trnsposed input table file to find total counts per k-mer: $!";
my @removePatterns;
my $totalPatterns=0;
$min = 100000;
$max = 0;

#put totals per k-mer into an array - this array is of size $totalPatterns
my @totalCountersPerpattern;
my $totalPatterns=0;
my $lineCounter=0;

print "\n".'_____________________________'."\n";
print "Started expected genomic copies calculations with Bayes model\n";
print '*****************************'."\n\n";
while (my $line = <$tableReader>)
{
   	chomp $line;
	my @lineArray = split(',', $line);
	my $total=0;

    #each line represents a variant
    $removePatterns[$lineCounter]=0; #initialize array entry for this variant
    for(my $i=0; $i <= $#lineArray; $i ++)
	{
        $total += abs ($lineArray[$i]);        
    }

	if($total > $max)
    {
        $max=$total;
    }
    if($total < $min)
    {
        $min=$total;
    }
	push (@totalCountersPerpattern, $total);
	
    
    if($totalIndividuals-1 != $#lineArray)
    {
        my $actual =  $#lineArray+1;
        print "Error in input tables: should be $totalIndividuals columns in the transposed table on line $lineCounter but there are in fact $actual columns \n";
        exit (1);
    }
    
    $lineCounter++;
}
close $tableReader;

$totalPatterns = $lineCounter;

print "Started with $totalIndividuals samples and $totalPatterns patterns\n";

print  "Total counts of k-mers for all samples range from  $min to $max \n";

#************************************
#3. For each k-mer i estimate total number of its genomic copies in the entire population
#here we also want to mark as invalid all k-mers with estimated copy number > 2*$totalIndividuals
#*******************************

my $lambdaDepthTotal = 0;
for (my $j=0; $j <= $#lambdaPerFile; $j ++) 
{
	$lambdaDepthTotal += $lambdaPerFile[$j];
}

my $meanDepth = $lambdaDepthTotal/(2*$totalIndividuals); 

print STDERR "Average sequencing depth= ".$meanDepth."\n";

my @totalCopiesPerPattern;
$min = 100000;
$max = 0;

for (my $i=0; $i < $totalPatterns; $i++) 
{
	my $estPatCopy = $totalCountersPerpattern[$i]/$meanDepth;
	if($estPatCopy > 2*$totalIndividuals)
	{
		$estPatCopy = 2*$totalIndividuals;
        $removePatterns[$i]=2; #2 stands for num of copies in all samples greater than expected for 2n
	}

	push (@totalCopiesPerPattern, $estPatCopy);
	if($estPatCopy > $max and !$removePatterns[$i])
	{
		$max=$estPatCopy;
	}
	if($estPatCopy < $min and !$removePatterns[$i])
	{
		$min=$estPatCopy;
	}
}
print "Estimated number of genetic copies per k-mer ranges from $min to $max "."\n";

#************************************
#4. For each k-mer i estimate allele frequency f 
#*******************************
my @alleleFrequencyPerPattern;

$min = 100000;
$max = 0;

for (my $i=0; $i < $totalPatterns; $i ++) 
{
	my $f = $totalCopiesPerPattern[$i]/(2*$totalIndividuals);
	push (@alleleFrequencyPerPattern,$f );
	if($f > $max and !$removePatterns[$i])
	{
		$max=$f;
	}
	if($f < $min and !$removePatterns[$i])
	{
		$min=$f;
	}
}

print "Allele frequency per k-mer (f) ranges from $min to $max "."\n";

#************************************
#5. for each k-mer i find probability of AA, A-, and -- (prior probabilities according to HW)
#*******************************
my $e = 2.71828;

my @priors;  #2d array - prior probabilities that this k-mer represents [0] homozygote AA, [1] heterozyg A-, or [2] homozyg --
$min = 1000000;
$max = 0;
my $log_2 = log 2;
for (my $i=0; $i < $totalPatterns; $i ++) 
{
    
    my $f = $alleleFrequencyPerPattern [$i];

    $priors [$i] [0] =(1-$f) *(1-$f); #P0 (1-f)^2
    $priors [$i] [1] =(2*$f)*(1-$f) ; #P1 2f(1-f)
    $priors [$i] [2] = $f*$f ; #P2 f^2	
    if($priors [$i] [0] > $max and !$removePatterns[$i] )
    {
	    $max=$priors [$i] [0];
    }
    if($priors [$i] [1] > $max and !$removePatterns[$i] )
    {
	    $max=$priors [$i] [1];
    }
    if($priors [$i] [2] > $max and !$removePatterns[$i]  )

    {
	    $max=$priors [$i] [2];
    }

    if($priors [$i] [0] < $min and !$removePatterns[$i] )
    {
	    $min=$priors [$i] [0];
    }
    if($priors [$i] [1] < $min and !$removePatterns[$i] )
    {
	    $min=$priors [$i] [1];
    }
    if($priors [$i] [2]< $min and !$removePatterns[$i] )
    {
	    $min=$priors [$i] [2];
    }
    
}

print "prior probabilities range from $min to $max "."\n";

#now we have @priors for each k-mer, @lambdaPerFile for each file, @totalCountersPerpattern, 
# calculating expected values, reading each row in the tableSamplesInRows first - 
# to 
#***********************************************************************
# 6. create a matrix file where each row is a sample, each column is an expected value of a k-mer's genomic copy
#***********************************************************************
open my $tableReader, $tableSamplesInRows or die "error trying to open input table file: $!";

my $minAlpha = 100000;
my $maxAlpha = 0;
$min = 1000000;
$max = 0;

#for each file and each pattern find conditional probabilities CP(gij=0|cij), log CP(gij=1|cij), log CP(gij=2|cij) and general probability P(cij)
#compute expected values and add them to the output line- then check and mark as 1 all columns where values fail
my $row=0; #sample id
while (my $line = <$tableReader>)
{
   	chomp $line;
	my @lineArray = split(',', $line);
	
	my $lambda= $lambdaPerFile[$row];
	my $halfLambda = $lambda/2;
	
	for (my $col=0; $col <= $#lineArray; $col ++) 
	{
		my $counter= abs (int ($lineArray[$col]));
		my @condProb;
		#compute CP0
		if($counter >0)
		{
			$condProb [0] = 0;
		}
		else
		{
			$condProb [0] = 1;
		}
		
		my $logFactorial = 0;
		for (my $i=2; $i <= $counter; $i ++)
		{
			$logFactorial+= log ($i); 
		}
		
		my $counterFactorial = $e ** $logFactorial;		

		#compute CP1
		$condProb  [1] = ($halfLambda ** $counter) * ($e ** (-$halfLambda)) /$counterFactorial;
		
		#compute CP2
		$condProb  [2] = ($lambda ** $counter) * ($e ** (-$lambda)) /$counterFactorial;		

		if($condProb  [1] > $max and !$removePatterns[$col] )
		{
			$max=$condProb  [1];
		}

		if($condProb  [1] < $min and !$removePatterns[$col] )
		{
			$min=$condProb  [1];
		}

		if($condProb [2] > $max and !$removePatterns[$col] )
		{
			$max=$condProb  [2];
		}

		if($condProb [2] < $min and !$removePatterns[$col] )
		{
			$min=$condProb  [2];
		}

		#compute alpha
		$condProb  [3] = $priors [$col] [0] * $condProb  [0] + $priors [$col] [1] * $condProb  [1] +$priors [$col] [2] * $condProb [2];		

		if($condProb  [3] < $minAlpha and !$removePatterns[$col] )
		{
			$minAlpha=$condProb  [3];
		}

		if($condProb  [3] > $maxAlpha and !$removePatterns[$col] )
		{
			$maxAlpha=$condProb  [3];
		}

		my $expect=0;
		if (!$removePatterns[$col])
		{
			if($condProb [3] != 0)
			{
				$expect = (2*$condProb  [2] *$priors [$col][2] + $condProb  [1]*$priors [$col][1])/$condProb [3];
				if($expect != $expect or $expect eq "inf" or !($expect < 9**9**9)) #check for invalid value
				{
					$removePatterns[$col]=3;
					$expect=0;
				}
			}
			else
			{
				$removePatterns[$col]=4;					
			}
		}
        if($col == 0) #first entry on this line
        {
            print TEMPWRITER $expect;
        }
        else
        {
		    print TEMPWRITER ",".$expect;
        }
	}
    $row++;
	#print "Processed row $row \n";
	print TEMPWRITER "\n";
	print "Written row $row to temp table \n";	
}
close $tableReader;
close TEMPWRITER;
print "Conditional probabilities range from $min to $max "."\n";
print "Alpha probabilities range from $minAlpha to $maxAlpha "."\n";



#*******************************
#9. Now we are going to perform the same procedure but now for the transposed table
#*****************************************
#here we can skip invalid rows (removepatterns true - right away)
open my $tableReader, $tableSamplesInColumns or die "error trying to open transposed input table file: $!";
#for each file (col) and each pattern (row) find conditional probabilities CP(gij=0|cij), log CP(gij=1|cij), log CP(gij=2|cij) and general probability P(cij)
#if a valid row (pattern) - compute expected values and add them to the output line-
my $row=0; #kmer id
my $validRow=0;
while (my $line = <$tableReader>)
{
   	if(!$removePatterns[$row])
    {
        chomp $line;
	    my @lineArray = split(',', $line);    
	
        my @expectedForThisKmer;
    
	    for (my $col=0; $col <= $#lineArray; $col ++) 
	    {
            my $lambda= $lambdaPerFile[$col];
	        my $halfLambda = $lambda/2;

		    my $counter= abs (int ($lineArray[$col]));
		    my @condProb;
		    #compute CP0
		    if($counter >0)
		    {
			    $condProb [0] = 0;
		    }
		    else
		    {
			    $condProb [0] = 1;
		    }
		
		    my $logFactorial = 0;
		    for (my $i=2; $i <= $counter; $i ++)
		    {
			    $logFactorial+= log ($i); 
		    }
		
		    my $counterFactorial = $e ** $logFactorial;		

		    #compute CP1
		    $condProb  [1] = ($halfLambda ** $counter) * ($e ** (-$halfLambda)) /$counterFactorial;
		
		    #compute CP2
		    $condProb  [2] = ($lambda ** $counter) * ($e ** (-$lambda)) /$counterFactorial;			

		    #compute alpha
		    $condProb  [3] = $priors [$row] [0] * $condProb  [0] + $priors [$row] [1] * $condProb  [1] +$priors [$row] [2] * $condProb [2];	
		
            my $expect=0;
		
			if($condProb [3] != 0)
			{
				$expect = (2*$condProb  [2] *$priors [$row][2] + $condProb  [1]*$priors [$row][1])/$condProb [3];
				if($expect != $expect or $expect eq "inf" or !($expect < 9**9**9)) #check for invalid value
				{
					print "Unexpected error: invalid count detected where only valid should remain\n";
					exit (1);
                }
			}
			else
			{
				print "Unexpected error: invalid count detected where only valid should remain\n";
				exit (1);				
			}

            push(@expectedForThisKmer, $expect);
           
		}

        #determine if variance is non-zero - not all values are the same across all the samples     
           
        my $zero_variance = 1;
        my $currExpected =int ($expectedForThisKmer[0]+0.5); #rounding to the closest integer - 0, 1 0r 2	 $expectedForThisKmer[0];
		
		for (my $r=1; $r < $totalIndividuals and $zero_variance; $r ++) 
		{
			if(int ($expectedForThisKmer[$r]+0.5) != $currExpected)
			{
				$zero_variance=0;
			}
			
		}
		if($zero_variance)
		{
			$removePatterns [$row] =5;
		}
        
        print WRITERBYCOL $expectedForThisKmer[0];
        for (my $r=1; $r < $totalIndividuals; $r ++) 
	    {
	        print WRITERBYCOL ",".$expectedForThisKmer[$r];
        }
        print WRITERBYCOL "\n";
        if($row % 1000 == 0)
        {
            print "Written variant $row to transposed table\n";
        }	
        $validRow++;
        
	}
    $row++;	
}
close $tableReader;
close WRITERBYCOL;

print "Finished calculating expected values for $validRow valid k-mers and $totalIndividuals samples in transposed table\n";

#we are iterating over final expected values and chop out all the columns which are marked as invalid
#***********************************************************************
# 8. Remove invalid columns from the original table 
#***********************************************************************
#now we are iterating over final expected values and chop out all the columns which are marked as invalid
$row=0;
open my $tableReader, "$tableSamplesInRows"."TEMP" or die "error trying to open temp output table file for reading: $!";
while (my $line = <$tableReader>)
{
   	chomp $line;
	my @lineArray = split(',', $line);

	#print "Extracted file name $fname \n";
	
	my $totalColumnsPerFile = 0;
	for ( my $col=0; $col <= $#lineArray; $col++) 
	{
		if (!$removePatterns[$col] or $removePatterns[$col]==5  )
		{
			my $expected= $lineArray[$col];
            if($totalColumnsPerFile == 0)
            {
			    print WRITERBYROW $expected;
            }
            else
            {
                print WRITERBYROW ','.$expected;
            }
			$totalColumnsPerFile++;
		}
	}
	print WRITERBYROW "\n";
	print "Remains $totalColumnsPerFile valid columns for file $row \n";
    $row++;
}

close $tableReader;
close WRITERBYROW;

print "\n".'*****************************'."\n";
print "Finished both tables with $validRow valid k-mers and $totalIndividuals samples\n";
print '_____________________________'."\n\n";

print REMOVEWRITER $removePatterns[0];
for (my $i=1; $i<=$#removePatterns; $i++)
{
    print REMOVEWRITER ','.$removePatterns[$i];
}
