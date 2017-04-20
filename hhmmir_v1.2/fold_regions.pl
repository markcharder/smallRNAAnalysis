#################################################
#
# Copyright @Sabah Kadri. Contact Sabah Kadri for use of this file outside of HHMMiR
#
#################################################


#!/usr/bin/perl

# execute as perl hairpin.pl <pathname_of_RNAfold> <fasta_file_of_genome>

# $ARGV[0] = pathname of RNAfold
# $ARGV[1] = fasta file with sequence data

# Reading command line arguments
$numArgs = $#ARGV + 1;
# If incorrect number of arguments, end program
if($numArgs < 2){
	die("Incorrect number of input parameters\n");
}

# Opening the fasta file
# The fourth argument will be the huge fa file
open(INFILE,$ARGV[1]) || die("Could not open file!");


# Reading data
#line at a time
my $ptr=0;
#my $window=300;
my $window=1000;
my $chrom_name;
my $kb_data="";
my $data;
chomp($data = <INFILE>);
my $inputfile = "ip.txt";
my $outputfile = "folded.fa";

# While there is more sequence data
while($data){
	
	# if the first line contains a > 
	# start reading scaffold 
	if(!((index($data,'>')==-1))){
		

		$chrom_name = $data;
		
		# open output file
		open (OUTPUTFILE, '>>'.$outputfile);
		# print chromosome name to output file	
		print OUTPUTFILE $chrom_name."\n";
		close(OUTPUTFILE);
		
		# The next line will naturally be the sequence
		chomp($data = <INFILE>);
		
		# till reach EOF or next scaffold
		while((index($data,'>')==-1) && ($data)){
			
			# read 1000 nts
			if((length($kb_data)+length($data))<$window){
				# if sequence contains no gaps i.e. N
				if(index($data,'N')==-1){
					$kb_data = $kb_data.$data;
					$ptr=0;
				}
				# if sequence contains gaps, terminate the sequence for hairpin fold there
				else{
					$ptr = index($data,'N');
					$kb_data = $kb_data.substr($data,0,$ptr);
					
					#if length of $kb_data is zero, dont fold
					if(length($kb_data)>0){
						
						# 1.Write it to input and output files
						open (INPUTFILE, '>'.$inputfile);
						print INPUTFILE $kb_data;
						print INPUTFILE "\n";
						close(INPUTFILE);
						
						# 2.Call the RNAfold program
							
						use strict;
						use warnings;
						# $ARGV[0] = pathname of RNAfold
						
						my $status = system("$ARGV[0] < $inputfile >> $outputfile");
						open (OUTPUTFILE, '>>'.$outputfile);	
						print OUTPUTFILE "\n";
						close(OUTPUTFILE);
					}
					
					$ptr = rindex($data,'N');
					
					if($ptr<(length($data)-1)){
						
						$kb_data="";
						$kb_data = $kb_data.substr($data,$ptr+1,length($data)-1-$ptr);
						#print(substr($data,$ptr+1,length($data)-1-$ptr));
						#print "\n";
						#print($ptr);
						#print "\n";
					}
					#if gaps extend till end of line, dont form a new $kb_data yet
					else{
						$kb_data="";
					}
					$ptr = 0;

				}
			}
			
			else{
				$ptr = $window - length($kb_data);
				$kb_data = $kb_data.substr($data,0,$ptr+1);
				#print(length($kb_data));
				#print $kb_data;
				#print "\n";
				# window size read
				# 1.Write it
				open (INPUTFILE, '>'.$inputfile);
				print INPUTFILE $kb_data;
				print INPUTFILE "\n";
				close(INPUTFILE);
				
				# 2.Call the RNAfold program
					
				use strict;
				use warnings;
				# $ARGV[0] = pathname of RNAfold
				# $ARGV[1] = inputfile name in curr. dir
				# $ARGV[2] = outputfile name in curr. dir
				my $status = system("$ARGV[0] < $inputfile >> $outputfile");
				open (OUTPUTFILE, '>>'.$outputfile);	
				print OUTPUTFILE "\n";
				close(OUTPUTFILE);
				
				$kb_data = substr($kb_data,(length($kb_data)-150)).substr($data,$ptr);
				$ptr=0;
	
			}
			chomp($data = <INFILE>);
			
		
		}
		
		
		# If the next line contains '>', it means we reached next scaffold
		# Or if we reach the EOF
		# Now perform the calculations on the kb_data
	
		if(!((index($data,'>')==-1)) || (!$data))
		{
			# 1.Write it
			open (INPUTFILE, '>'.$inputfile);
			print INPUTFILE $kb_data;
			print INPUTFILE "\n";
			close(INPUTFILE);
			
			# 2.Call the RNAfold program
				
			use strict;
			use warnings;
			# $ARGV[0] = pathname of RNAfold
			# $ARGV[1] = inputfile name in curr. dir
			# $ARGV[2] = outputfile name in curr. dir
			my $status = system("$ARGV[0] < $inputfile >> $outputfile");
			open (OUTPUTFILE, '>>'.$outputfile);	
			print OUTPUTFILE "\n";
			close(OUTPUTFILE);

			$kb_data = "";
			$ptr=0;
			
		}
		
	}
	else{
		chomp($data = <INFILE>);
	}
	
	# read 1KB into a string
	# or till end of scaffold

}


# Close File
close(INFILE); 
