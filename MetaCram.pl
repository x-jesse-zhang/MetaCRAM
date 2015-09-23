#!/usr/bin/perl

# use lib '~/perl5/lib/perl5';

use strict;
use warnings;
use File::Find;
use File::Slurp;
use File::stat;
use File::Copy;
use Array::Utils qw(:all);
use Getopt::Long;
use File::Basename;
use Capture::Tiny ':all';
use Cwd;
# use Statistics::R;
use Config;
use FileHandle;
use threads;
use Thread::Semaphore;
use Data::Dumper;
use Cwd;
use Tie::File;

STDOUT->autoflush(1);

my $PROG = basename $0;

$Config{useithreads} or die "Your compilation of perl does not support threads.\n";

# Help stuff
sub usage {
	my $exit_code = @_ ? shift : 64;
	print STDERR <<EOF;
Usage: $PROG [options] <filename(s)>

Options:
   Compression (--compress) :
       --db NAME            Absolute path to kraken's database
       --blast_db_location  Absolute path to location of BLAST's database
       --threads NUM        Number of threads to use
       --output NAME        Name for MCUIUC's output folder
       --kraken_dir NAME    Absolute path to the folder containing kraken, not kraken itself
       --bowtie_dir NAME    Absolute path to the folder containing bowtie2, not bowtie2 itself
       --blast_dir NAME     Absolute path to the folder containing BLAST, namely blastn
       --genus_db NAME      Absolute path to where NCBI genus database is
       --paired             The two filenames provided are paired-end reads  
       --cram               Absolute path to Cramtools' jar file
       --preserve-all       Preserve Bowtie's alignment rate, species identified by Kraken and BLAST, 
                            contig files from IDBA_UD and SAM files from Bowtie
    
       --preserve-alignment-rate     Preserve Bowtie's alignment rate
       --preserve_species   Preserve species identified by Kraken and BLAST
       --preserve-contig	Preserve contig files from IDBA
       --preserve-sam       Preserve SAM files from Bowtie
       --huffman            Use Huffman encoding for CRAM
       --golomb             Use Golomb encoding for CRAM
       --exGolomb           Use Extended Golomb encoding for CRAM

   Decompression (--decompress) :
   		Work in progress.
EOF
	exit $exit_code;
}

sub display_help {
	usage(0);
}

#Settings
my $compress;
my $decompress;
#dev: default 4
my $threads;
my $paired='';
#dev: McuiucOutput
my $output_dir;
#dev: default '/ramdisk/minikraken_20141208'
my $kraken_db_location;
#dev: default '/opt/kraken'
my $kraken_dir;
#dev: default '/shared/genus_db2'
my $genus_db_location;
#dev: default '/opt/bowtie2'
my $bowtie2_location;
#dev: default '/opt/ncbi-blast-2.2.29+/bin'
my $blast_location;
#dev: default '/opt/ncbi-blast-2.2.29+/db'
my $blast_db_location;
#dev: default '/shared3/cramtools-2.1.jar'
# my $cram;
#dev: default 0
my $preserve_bowtie = 0;
my $preserve_species = 0;
my $preserve_contig = 0;
my $preserve_sam = 0;
my $preserve_all;

my $huffman = 0;
my $golomb = 0;
my $exGolomb = 0;
my $huffmanPath = getcwd . "/HuffmanCRAM.jar";
my $golombPath = getcwd . "/GolombCRAM.jar";
my $exGolombPath = getcwd . "/ExtendedGolombCRAM.jar";
my $mfcPath = getcwd . "/MFCompressC";
my $mfdPath = getcwd . "/MFCompressD";
my $maxRounds = 2;

my $fna;
my $input1;
my $input2;
my $totalReads;
my @SAMlines;
my $alignFile;
my @SAMfilteredFiles;
my @genus;
my $fileID;
my @unaligned_parts;
my $idba_splited = 0;
my $roundOfUnaligned = 1;
my $splitCount;
my $bowtieOut;
my $SAMfile;
my @headers;
my $readMap_dir;
my %startposHash;
my %fnaHash;
my $semaphore = Thread::Semaphore->new();

# TODO, modify this so it's optional to spcify directories and rely on native PATH
# Command line options
GetOptions(
	"compress" => \$compress,
	"decompress" => \$decompress,
	"help" => \&display_help,
	"db=s" => \$kraken_db_location,
	"threads=i" => \$threads,
	"output=s" => \$output_dir,
	"kraken_dir=s" => \$kraken_dir,
	"paired" => \$paired,
	"genus_db=s" => \$genus_db_location,
	"bowtie_dir=s" => \$bowtie2_location,
	"blast_dir=s" => \$blast_location,
	"blast_db=s" => \$blast_db_location,
	# "cram=s" => \$cram,
	"preserve-alignment-rate" => \$preserve_bowtie,
	"preserve-contig" => \$preserve_contig,
	"preserve-sam" => \$preserve_sam,
	"preserve-all" => \$preserve_all,
	"huffman" => \$huffman,
	"golomb" => \$golomb,
	"exGolomb" => \$exGolomb
);
if ($compress && $decompress) {
	print "Cannot specify both compress and decompress.\n";
	exit(0);
}

if (! defined $compress && ! defined $decompress) {
	usage(0);
}

if ($decompress) {
	do_decompress();
}

#If conditionals to gather the necessary information from terminal
if (! defined $threads) {
	#print("Default to using 4 threads.\n");
	$threads=4;
}

if (! defined $output_dir) {
	#print("Default to output to MetaCRAM\n");
	$output_dir="MetaCRAM";
}

if (! defined $kraken_db_location) {
	#die "must specify database name with --db\n";
	$kraken_db_location="/ramdisk/minikraken_20140330";
}

if (! defined $kraken_dir) {
	#die "must specify kraken directory\n";
	$kraken_dir='/opt/kraken';
}

if (! defined $genus_db_location) {
	#die "must specify genus database location\n";
	$genus_db_location='/shared/genus_db2';
}

if (! defined $bowtie2_location) {
	#die "must specify bowtie2 directory\n";
	$bowtie2_location="/opt/bowtie2";
}

if (! defined $blast_location) {
	#die "must specify blast directory\n";
	$blast_location="/opt/ncbi-blast-2.2.29+/bin";
}

if (! defined $blast_db_location) {
	#die "must specify blast database location\n";
	$blast_db_location="/opt/ncbi-blast-2.2.29+/db";
}

# if (! defined $cram) {
# 	#die "must specify cramtools location\n";
# 	$cram="/shared3/ExtendedGolombCRAM.jar";
# }

if ($preserve_all) {
	$preserve_contig = 1;
	$preserve_species = 1;
	$preserve_bowtie = 1;
	$preserve_sam = 1;
}

if ($huffman == 0 && $golomb == 0 && $exGolomb == 0) {
	die "Must specify a compression method\n";
}

if ($huffman + $golomb + $exGolomb > 1) {
	die "Cannot specify more than 1 compression method\n";
}

if (! -e $huffmanPath) {
	die "Cannot find HuffmanCRAM.jar in running directory of MetaCRAM.pl\n";
}

if (! -e $golombPath) {
	die "Cannot find GolombCRAM.jar in running directory of MetaCRAM.pl\n";
}

if (! -e $exGolombPath) {
	die "Cannot find ExtendedGolombCRAM.jar in running directory of MetaCRAM.pl\n";
}

if (! @ARGV) {
	print STDERR "Need to specify input files\n";
	usage();
}

if ($paired && @ARGV != 2) {
	#die "--paired requires exactly two filenames!";
	my ($filename, $directories, $suffix) = fileparse("$ARGV[0]", qr/\.[^.]*/);
	my $firstPart = (split(/\./, $ARGV[0]))[0];
	$input1= $firstPart . "_1" . $suffix;
	$input2= $firstPart . "_2" . $suffix;
	# print("--paired but only 1 input file.\n");		

	# unless (-e 'sepReads.class') {
	# 	print("Compiling sepReads.java. ");
	# 	system("javac sepReads.java");
	# 	print("Done\n");
	# }

	print("Splitting $ARGV[0] into $input1 and $input2. ");
	# system("java sepReads $ARGV[0] $input1 $input2");
	sepReads($ARGV[0], $input1, $input2);
	print("Done.\n");
} elsif ($paired && @ARGV == 2) {
	$input1 = $ARGV[0];
	$input2 = $ARGV[1];
}

sub make_blast_config {
	open(my $config, ">", ".ncbirc");
	print $config "; Start the section for BLAST configuration\n";
	print $config "[BLAST]\n";
	print $config "; Specifies the path where BLAST databases are installed\n";
	print $config "BLASTDB=$blast_db_location\n";
	close $config;
}

# print("Counting amount of reads in input. ");
# if ($paired) {
# 	open(INPUT, $input1) or die "Can't open $input1";
# 	while(<INPUT>) {
# 		$totalReads++;
# 	}
# 	$totalReads = $totalReads / 2;
# 	close INPUT;
# } else {
# 	open(INPUT, $ARGV[0]) or die "Can't open $ARGV[0]";
# 	while(<INPUT>) {
# 		$totalReads++;
# 	}
# 	close INPUT;
# 	$totalReads = $totalReads / 4;
# }
# print("Done\n");


# More definition
my $output_dir_temp="$output_dir/MetaCRAM";
my $reportFiltered="$output_dir_temp/reportFiltered";

my $kraken="$kraken_dir/kraken";
my $kraken_mpa_report="$kraken_dir/kraken-mpa-report";
my $kraken_output_dir="$output_dir_temp/krakenOutput";

my $bowtie_output_dir="$output_dir_temp/bowtieOutput";

my $logFile="$output_dir/McuiucLog";

my $SAMfilter_dir="$output_dir_temp/SAMFilter";

my $idbaAssembled_dir="$output_dir_temp/assembled";

my $blast_output_dir="$output_dir_temp/blastOUT";

my $cram_dir="$output_dir_temp/cramOutput";

my $min_contig_length=10**6; # Length of the shortest contig to which to consider for BLAST identification
my $twentyGB=21474836480;
my $fiveGB=5300000000;
#my $max_assembly_size=(10**8); # 95 MB
my $max_assembly_size=2147483648; # 2GB
my $home_dir = getcwd();

# COMMAND LINE INPUT DEBUGGING
# print("\n");
# print("$kraken_db_location" . "\n");
# print("$threads" . "\n");
# print("$output_dir" . "\n");
# print("$kraken_dir" . "\n");
# print("$paired" . "\n");
# print("$genus_db_location" . "\n");
# print("$bowtie2_location" . "\n");
# if ($paired) {
# 	print($ARGV[0] . " " . $ARGV[1] . "\n\n");
# } else {
# 	print($ARGV[0] . "\n\n");
# }

# Making the necessary folders
print("Making necessary directories. ");
unless (-d $output_dir) {
	mkdir($output_dir) or die "Output directory cannot be created.\n";
}

unless (-d $output_dir_temp) {
	mkdir($output_dir_temp) or die "Output directory cannot be created.\n";
}

unless (-d $SAMfilter_dir) {
	mkdir($SAMfilter_dir) or die "Output directory cannot be created.\n";
}

unless (-d $kraken_output_dir) {
	mkdir($kraken_output_dir) or die "Output directory cannot be created.\n";
}

unless (-d $bowtie_output_dir) {
	mkdir($bowtie_output_dir) or die "Output directory cannot be created.\n";
}

unless (-d $SAMfilter_dir) {
	mkdir($SAMfilter_dir) or die "Output directory cannot be created.\n";
}

# unless (-d $readMap_dir) {
	# mkdir($readMap_dir) or die "Output directory cannot be created.";
# }

unless (-d $idbaAssembled_dir) {
	mkdir($idbaAssembled_dir) or die "Output directory cannot be created.\n";
}

unless (-d $blast_output_dir) {
	mkdir($blast_output_dir) or die "Output directory cannot be created.\n";
}

unless (-d $cram_dir) {
	mkdir($cram_dir) or die "Output directory cannot be created.\n";
}

unless (-e ".ncbirc") {
	make_blast_config();
}

print("Done.\n");

# unless (-e 'sepReads.class') {
# 	print("Compiling sepReads.java");
# 	system("javac sepReads.java");
# 	print "Done.\n"
# }

# Adding tags to second mate file
# print "Tagging the second mate file. ";
# my $originalInput2 = $input2;
# {
# 	my @parts = split("/", $input2);
# 	my $name = (split(/\./, $parts[-1]))[0];
# 	my $extension = (split(/\./, $parts[-1]))[1];
# 	my $path = join("/", @parts[0..(@parts - 2)]);

# 	$input2 = "$path/modified$name.$extension";

# 	open(my $input, "<", $originalInput2) or die "Cannot open $originalInput2\n";
# 	open(my $output, ">", $input2) or die "Cannot create a modified version of $input2\n";

# 	while (my $line = <$input>) {
# 		chomp $line;

# 		if (substr($line, 0, 1) ne ">") {
# 			print $output $line . "\n";
# 			next;
# 		}

# 		@parts = split(" ", $line);
# 		$parts[0] = $parts[0] . '\2';

# 		$line = join(" ", @parts);
# 		print $output $line . "\n";
# 	}

# 	close $input;
# 	close $output;
# }
# print "Done.\n";

sub run_kraken {
	if (-e "$kraken_output_dir/.cpKraken") {
		print "Kraken checkpoint file exists. Skipping.\n";
		return;
	}

	print("Running Kraken and kraken_mpa_report. ");

	if ($paired) {
		system("$kraken --db $kraken_db_location --threads $threads --paired $input1 $input2>$kraken_output_dir/KrakenOutput");
	} else {
		system("$kraken --db $kraken_db_location --threads $threads $ARGV[0]>$kraken_output_dir/KrakenOutput");
	}

	system("$kraken_mpa_report --db $kraken_db_location $kraken_output_dir/KrakenOutput>$kraken_output_dir/KrakenOutputReport");

	open(my $cp_kraken, ">", "$kraken_output_dir/.cpKraken") or die "Kraken checkpoint file cannot be created.\n";
	close $cp_kraken;

	print("Done.\n");
}

sub filter_kraken {
	my %fileSizes;
	my @fileSorted;
	my @f3;
	my @specieNames;
	my @fnaFiles;
	my $filteredReads;
	my @mpaReportLines;
	my %species;
	my @speciesSorted;
	my $fnaForPrint;

	if (-e "$kraken_output_dir/listOfFnaFiles") {
		print "Kraken filter checkpoint file exists. Skipping.\n";
		my @fnas = read_file("$kraken_output_dir/listOfFnaFiles", chomp => 1);
		$fna = join(",", @fnas);
		return;
	}

	print "Filtering Kraken output.\n";

	@mpaReportLines = read_file("$kraken_output_dir/KrakenOutputReport", chomp => 1);

	foreach (@mpaReportLines) {
		my ($specie, $reads) = split("\t");
		# $totalReads += $reads;
		my @fullName = split(/\|/, $specie);

		if (! $fullName[6]) {
			next;
		}

		$filteredReads += $reads;

		$species{substr($fullName[6], 3)} = $reads;
	}

	@speciesSorted = sort {$species{$b} <=> $species{$a}} keys %species;

	# static threshold = bad
	if (scalar(@speciesSorted) > 75) {
		@specieNames = @speciesSorted[0..74];
	} else {
		@specieNames = @speciesSorted;
	}

	# Getting .fna files for bowtie
	print("Gathering .fna files\n");

	#@specieNames = read_file($reportFiltered, chomp => 1);

	foreach my $name (@specieNames) {
		#print("$name\n");
		my @folder_names = `ls $genus_db_location | grep $name`;
		#my $alternative_name;

		# if specie is not found, split name into parts and look up using the possible 3rd place variation name
		if (! @folder_names) {
			print("Using variation to find $name\n");
			my @name_parts=split('_', $name);

			@folder_names = `ls $genus_db_location | grep $name_parts[-1]`;

			# if still not found, look up the genus
			if (! @folder_names) {
				print("Using genus to find $name\n");
				@folder_names = `ls $genus_db_location | grep $name_parts[0]`;

				# if genus doesn't exist, ask user
				if (! @folder_names) {
					print("\nCan't find $name in genus database.\nIs there an alternative genus? ");
					$name=<STDIN>;
					chomp($name);
					print("\n");
				}
			}

			chomp(@folder_names);
			#print join("\n", @folder_names), "\n";


			# print("\nCan't find $name in genus database.\nIs there an alternative genus? ");
			# #$alternative_name=<STDIN>;
			# $name=<STDIN>;
			# chomp($name);

			# #@folder_names = `ls $genus_db_location | grep $alternative_name`;
			# @folder_names = `ls $genus_db_location | grep $name`;
		}

		my @specieFiles;

	# ERROR: if folder contains subfolders, that subfolder will get listed as an .fna file
		foreach my $folder (@folder_names) {
			chomp($folder);
			my @files = read_dir("$genus_db_location/$folder", prefix => 1);
			chomp(@files);

			for (my $i = 0; $i < (scalar(@files)); $i++) {
				my @splitName = split("/", $files[$i]);
				@splitName = split(/\./, $splitName[@splitName - 1]);

				if (! $splitName[1] || $splitName[1] ne 'fna') {
					splice(@files, $i, 1);				
				}
			}

			push(@specieFiles, @files);

			# my %fileSizes = map {$_, -s} @files;
			# my @fileSorted = sort {$fileSizes{$b} <=> $fileSizes{$a}} keys %fileSizes;
			# my @f3;

			# if (@fileSorted < 3) {
			# 	@f3 = @fileSorted[0..(@fileSorted-1)];
			# } else {
			# 	@f3 = @fileSorted[0..2];
			# }

			#print (scalar @f3 . "\n");

			#my @ftable = map {{$_, -s}} @files;
			#my @fts = sort {$$b[1]-$$a[1]} @ftable;
			#my @ft3 = @fts[0..2];

			#push(@fnaFiles, @f3);
			#print join("\n", @f3), "\n";
			#<STDIN>;
			#push(@fnaFiles, @files);
			#print join("\n", @files), "\n";
		}

		%fileSizes = map {$_, -s} @specieFiles;
		@fileSorted = sort {$fileSizes{$b} <=> $fileSizes{$a}} keys %fileSizes;

		my $staticLimit = 5;

		if (@fileSorted < $staticLimit) {
			@f3 = @fileSorted[0..(@fileSorted-1)];
		} else {
			@f3 = @fileSorted[0..($staticLimit - 1)];
		}

		# print ("Total: " . scalar @specieFiles . "\n");

		# foreach (@fileSorted) {
		# 	printf "$_ %i\n", -s;
		# }

		# print "\nCollected: " . scalar @f3 . "\n";

		# foreach (@f3) {
		# 	printf "$_ %i\n", -s;
		# }

		push(@fnaFiles, @f3);
		#print join("\n", @f3), "\n";
		#<STDIN>;

		# print("\n");

		#chomp(@folder_names);

		#print($name . "\n");
		#print join("\n", @folder_names), "\n\n";
	} 

	print("Filter complete.\n");

	$fna=join(",", @fnaFiles);
	$fnaForPrint=join("\n", @fnaFiles);
	#print($fna);

	open(my $out, ">", "$kraken_output_dir/listOfFnaFiles");
	print $out $fnaForPrint;
	close $out;
}

# TODO BT2log is not working
sub run_bowtie {
	my $BTbuildLog;
	my $BT2Log;
	my $bowtieInput1;
	my $bowtieInput2;

	unless (-d "$bowtie_output_dir/round$roundOfUnaligned") {
		mkdir("$bowtie_output_dir/round$roundOfUnaligned") or die "Bowtie round$roundOfUnaligned output directory cannot be created.\n";
	}

	$bowtieOut="$bowtie_output_dir/round$roundOfUnaligned/BTout";
	$SAMfile="$bowtie_output_dir/round$roundOfUnaligned/BTSAMout";
	$BTbuildLog="$bowtie_output_dir/round$roundOfUnaligned/BTbuildLog";
	$BT2Log="$bowtie_output_dir/round$roundOfUnaligned/BT2Log";
	$SAMfilter_dir = "$output_dir_temp/SAMFilter/round$roundOfUnaligned";
	$readMap_dir = "$SAMfilter_dir/readMap";
	$cram_dir = "$output_dir_temp/cramOutput/round$roundOfUnaligned";
	%fnaHash = ();

	if ($roundOfUnaligned != 1) {
		# mkdir("$bowtie_output_dir/round$roundOfUnaligned") or die "Bowtie output directory cannot be created.";
		# $bowtieOut="$bowtie_output_dir/round$roundOfUnaligned/BTout";
		# $SAMfile="$bowtie_output_dir/round$roundOfUnaligned/BTSAMout";
		# $BTbuildLog="$bowtie_output_dir/round$roundOfUnaligned/BTbuildLog";
		# $BT2Log="$bowtie_output_dir/round$roundOfUnaligned/BT2Log";
		# if ($idba_splited) {
		# 	for (my $i = 1; $i <= $splitCount; $i++) {
		# 		$fna = join(",", $fna, "$blast_output_dir/round$roundOfUnaligned/part$i/out.fna");
		# 	}
		# } else {
			$fna="$blast_output_dir/round" . ($roundOfUnaligned - 1) . "/out.fna";
		# }
		$bowtieInput1 = "$output_dir_temp/SAMFilter/round" . ($roundOfUnaligned - 1) . "/unaligned" . ($roundOfUnaligned - 1) . "_1.fna";
		$bowtieInput2 = "$output_dir_temp/SAMFilter/round" . ($roundOfUnaligned - 1) . "/unaligned" . ($roundOfUnaligned - 1) . "_2.fna";
	} else {
		# mkdir("$bowtie_output_dir/initial") or die "Bowtie output directory cannot be created.";
		# $bowtieOut="$bowtie_output_dir/initial/BTout";
		# $SAMfile="$bowtie_output_dir/initial/BTSAMout";
		# $BTbuildLog="$bowtie_output_dir/initial/BTbuildLog";
		# $BT2Log="$bowtie_output_dir/initial/BT2log";
		$bowtieInput1 = $input1;
		$bowtieInput2 = $input2;
	}

	if (-e "$bowtie_output_dir/round$roundOfUnaligned/.cpBowtie") {
		print "Bowtie round$roundOfUnaligned checkpoint file exists. Skipping.\n";
		return;
	}

	print ("Running bowtie2-build. ");
	# DEBUG: disable this
	system("$bowtie2_location/bowtie2-build -f $fna $bowtieOut>$BTbuildLog");
	print("Done.\n");

	print("Running bowtie2. ");
	#REDIRECT OUTPUT TO STDOUT
	# DEBUG: disable this
	system("$bowtie2_location/bowtie2 --reorder --threads $threads --mm -x $bowtieOut -1 $bowtieInput1 -2 $bowtieInput2 -S $SAMfile -f --no-hd --no-sq>$BT2Log");
	print("Done.\n");

	open(my $cpBowtie, ">", "$bowtie_output_dir/round$roundOfUnaligned/.cpBowtie") or die "Bowtie round$roundOfUnaligned checkpoint file cannot be created.\n";
	close $cpBowtie;

	separateSAM();
}

sub read_start_position {
	my $pairNumber;
	my $debugSAMline;
	my @entry;
	my @previousEntry;
	my $notUnaligned;

	# print("SAMfilter dir is: $SAMfilter_dir\n");
	# print("SamFile is: $SAMfile\n");
	# <STDIN>;

	unless (-d $SAMfilter_dir) {
		mkdir($SAMfilter_dir) or die "Output directory cannot be created.\n";
	}

	if (-e "$SAMfilter_dir/.cpRSP") {
		print "Read start position round$roundOfUnaligned checkpoint file exists. Skipping.\n";
		return;
	}

	print("Reading start positions. ");

	# Replacement for filter.awk
	#@SAMlines=read_file($SAMfile, chomp => 1);

	print "Opening SAM file. ";
	open(SAMhandle, '<', $SAMfile) or die $!;
	print "Done.\n";

	print "Making unaligned file. ";
	open(my $unaligned, ">", "$SAMfilter_dir/unaligned$roundOfUnaligned") or die $!;
	print "Done.\n";

	print "Making unpaired file. ";
	open(my $unpaired, ">", "$SAMfilter_dir/unpaired") or die $!;
	print "Done.\n";

	print "Filtering SAM file. ";

	while (<SAMhandle>) {
		$notUnaligned = 0;

		if (@entry) {
			@previousEntry = @entry;
		}
		# } else {
		# 	$previousEntry[0] = "impossible-999";
		# }

		#print "In the loop\n";
		$debugSAMline++;
		# print;
		# print "\n";
		chomp;
		@entry = split("\t");
		# print join("\n", @entry), "\n\n";
		# <STDIN>;

		if (substr($entry[0], 0, 1) eq "@") {
			next;
		}

		# 1 flag not present means the read is unpaired
		if (($entry[1] & 1) == 0) {
			print $unpaired ">$entry[0]\n";
			print $unpaired "$entry[9]\n";
			next;
		}

		# 4 flag means "The read has no reported alignments"
		if (($entry[1] & 4) && $entry[2] eq '*') {
			$pairNumber = 0;

			if ($entry[1] & 128) {
				$pairNumber = 1;
			}

			if ($roundOfUnaligned == 1) {
				print $unaligned ">$entry[0]/$pairNumber\n";
			} else {
				print $unaligned ">$entry[0]\n";
			}
			
			print $unaligned "$entry[9]\n";
			next;
		}

		# my $tempName1;
		# my $tempName2;

		# if (substr($entry[0], -2, -1) eq "_" || substr($entry[0], -2, -1) eq "/") {
		# 	$tempName1 = substr($entry[0], 0, -2);
		# } else {
		# 	$tempName1 = $entry[0];
		# }

		# if (substr($previousEntry[0], -2, -1) eq "_" || substr($previousEntry[0], -2, -1) eq "/") {
		# 	$tempName2 = substr($previousEntry[0], 0, -2);
		# } else {
		# 	$tempName2 = $previousEntry[0];
		# }

		# # if this part runs, means read is aligned and paired
		# # if ($entry[0] eq $previousEntry[0]) {
		# if ($tempName1 eq $tempName2) {
		# 	# both are aligned
		# 	if ((($entry[1] & 4) == 0) && (($previousEntry[1] & 4) == 0)) {
		# 		# print("Both are aligned.\n");
		# 		# system("clear");
		# 		# print join("\n", @entry), "\n\n";
		# 		# print join("\n", @previousEntry), "\n\n";
		# 		# <STDIN>;
		# 		my ($AS, $useless, $score) = split(/:/, $entry[11]);
		# 		my ($oldAS, $oldUseless, $oldScore) = split(/:/, $previousEntry[11]);

		# 		if ($oldScore > $score) {
		# 			$alignFile = (split(/\|/, $previousEntry[2]))[3];
		# 		} else {
		# 			$alignFile = (split(/\|/, $entry[2]))[3];
		# 		}
		# 	} elsif ($previousEntry[1] & 4) { # if only the current mate is aligned
		# 		# print("Current m8 is aligned.\n");
		# 		$alignFile = (split(/\|/, $entry[2]))[3];
		# 		# print("entry[2] is $entry[2]\n");
		# 	} elsif ($entry[1] & 4) { #if only the other mate is aligned
		# 		# print("Other m8 is aligned.\n");
		# 		$alignFile = (split(/\|/, $previousEntry[2]))[3];
		# 		# print("previousEntry[2] is $previousEntry[2]\n");
		# 	} else { #if neither is aligned
		# 		die "Both gave alignment reference for both mates but there was no alignment.";
		# 	}

		# 	# print("New alignFile: $alignFile\n");
		# 	# <STDIN>;

		# 	$alignFile = (split(/\./, $alignFile))[0];
		# 	# print("New alignFile: $alignFile\n");
		# 	# <STDIN>;

		# 	open(my $alignments, ">>", "$SAMfilter_dir/$alignFile.alignments");
		# 	open(my $startpos, ">>", "$SAMfilter_dir/$alignFile.startpos");

		# 	if (! $startposHash{"$SAMfilter_dir/$alignFile.startpos"}) {
		# 		$startposHash{"$SAMfilter_dir/$alignFile.startpos"}++;
		# 	}

		# 	if ($entry[1] & 64) { # if current entry is first mate
		# 		print $alignments ">$entry[0]/0\n";
		# 		print $alignments "$entry[9]\n";
		# 		print $startpos "$entry[3]\n";

		# 		print $alignments ">$previousEntry[0]/1\n";
		# 		print $alignments "$previousEntry[9]\n";
		
		# 		print $startpos "$previousEntry[3]\n";
		# 	} else {
		# 		print $alignments ">$previousEntry[0]/0\n";
		# 		print $alignments "$previousEntry[9]\n";
		# 		print $startpos "$previousEntry[3]\n";

		# 		print $alignments ">$entry[0]/1\n";
		# 		print $alignments "$entry[9]\n";

		# 		print $startpos "$entry[3]\n";
		# 	}

		# 	close $alignments;
		# 	close $startpos;
		# }
	}

	close $SAMfile;
	close $unaligned;
	close $unpaired;

	print "Done.\n";

	if ($paired) {
		chdir($home_dir);
		print "Splitting unaligned into 2 files. ";
		sepReads("$SAMfilter_dir/unaligned$roundOfUnaligned", "$SAMfilter_dir/unaligned" . $roundOfUnaligned . "_1.fna", "$SAMfilter_dir/unaligned" . $roundOfUnaligned . "_2.fna");
		print "Done.\n";
	}

	open(my $cpRSP, ">", "$SAMfilter_dir/.cpRSP") or die "Read start position checkpoint file cannot be created.\n";
	close $cpRSP;
}

sub sepReads {
	my $mateSwitch = 0;
	my ($in, $out1, $out2) = @_;

	# open(INPUT, "$SAMfilter_dir/unaligned") or die "Cannot open unaligned file.";
	# open(my $output1, ">", "$SAMfilter_dir/unaligned_1.fna") or die "Cannot create output file.";
	# open(my $output2, ">", "$SAMfilter_dir/unaligned_2.fna") or die "Cannot create output file.";

	open(INPUT, $in) or die "Cannot open unaligned file.\n";
	open(my $output1, ">", $out1) or die "Cannot create output file.\n";
	open(my $output2, ">", $out2) or die "Cannot create output file.\n";

	while (my $line = <INPUT>) {
		chomp $line;

		if (! $mateSwitch) {
			print $output1 $line . "\n";
			$line = <INPUT>;
			chomp $line;
			print $output1 $line . "\n";
			$mateSwitch = 1;
		} else {
			print $output2 $line . "\n";
			$line = <INPUT>;
			chomp $line;
			print $output2 $line . "\n";
			$mateSwitch = 0;
		}
	}

	close INPUT;
	close $output1;
	close $output2;
}

sub getHeaders {
	print "Getting headers from output of blastdbcmd. ";
	open (my $fh, "<", "$blast_output_dir/round" . ($roundOfUnaligned - 1) . "/out.fna") or die "Can't open blastdbcmd output\n";

	for my $line (<$fh>) {
		chomp $line;
		# print $line . "\n";

		if (substr($line, 0, 1) ne ">") {
			next;
		}

		# <STDIN>;

		$line = substr($line, 1);

		# print "Pushing entry to headers\n";
		push(@headers, $line);
	}

	close $fh;

	print "Done.\n";
}

sub run_IDBA {
	chdir $SAMfilter_dir;
	$idba_splited = 0;

	if (-e "$idbaAssembled_dir/round$roundOfUnaligned/.cpIDBA") {
		print "IDBA checkpoint file exists. Skipping.\n";
		return;
	}

	print "Running IDBA_UD. ";

	# DEBUG: disable this
	unless (-d "$idbaAssembled_dir/round$roundOfUnaligned") {
		mkdir("$idbaAssembled_dir/round$roundOfUnaligned") or die "Output directory cannot be created.\n";
	}

	# TODO multithread this
	# if (-s "unaligned" >= $fiveGB) {
	# 	print "Unaligned file is too big. Splitting.\n";
	# 	system("split -b $max_assembly_size unaligned unaligned_parts");
	# 	@unaligned_parts=<$SAMfilter_dir/unaligned_parts*>;
	# 	$idba_splited = 1;

	# 	for my $i (@unaligned_parts) {
	# 		# DEBUG: disable this
	# 		mkdir("$idbaAssembled_dir/round$roundOfUnaligned/\$(basename $i)");
	# 		system("idba_ud --num_threads $threads -r $i -o $idbaAssembled_dir/round$roundOfUnaligned/\$(basename $i)");
	# 	}
	# } else {
		# DEBUG: disable this
		system("idba_ud --num_threads $threads -r $SAMfilter_dir/unaligned$roundOfUnaligned -o $idbaAssembled_dir/round$roundOfUnaligned/ > $idbaAssembled_dir/round$roundOfUnaligned/IDBALog");
	# }

	print "Done.\n";

	open(my $cp, ">", "$idbaAssembled_dir/round$roundOfUnaligned/.cpIDBA") or die "IDBA checkpoint file cannot be created.\n";
	close $cp;
}

sub run_BLAST {
	chdir $home_dir;

	if (-e "$blast_output_dir/round$roundOfUnaligned/.cpBLAST") {
		print "BLAST checkpoint file exists. Skipping.";
		return;
	}

	print "BLASTing.\n";

	# specify blast out name
	$splitCount = 0;

	unless (-d "$blast_output_dir/round$roundOfUnaligned") {
		mkdir("$blast_output_dir/round$roundOfUnaligned") or die $!;
	}

	# if ($idba_splited) {
	# 	for my $i (@unaligned_parts) {
	# 		$splitCount++;
	# 		# DEBUG: disable this
	# 		mkdir("$blast_output_dir/round$roundOfUnaligned/part$splitCount") or die "Output directory cannot be created.\n";
	# 		system("$blast_location/blastn -db nt -query $idbaAssembled_dir/round$roundOfUnaligned/\$(basename $i)/contig.fa -num_threads $threads -out $blast_output_dir/round$roundOfUnaligned/part$splitCount -max_target_seqs 1 -outfmt \"6 qseqid sseqid sgi sacc stitle evalue\"");
	# 	}
	# } else {
		# DEBUG: disable this
		system("$blast_location/blastn -db nt -query $output_dir_temp/assembled/round$roundOfUnaligned/contig.fa -num_threads $threads -out $blast_output_dir/round$roundOfUnaligned/blastOUT -max_target_seqs 1 -outfmt \"6 qseqid sseqid sgi sacc stitle evalue\"");
	# }

	print "Done.\n";

	open(my $cp, ">", "$blast_output_dir/round$roundOfUnaligned/.cpBLAST") or die "BLAST checkpoint file cannot be created.\n";
	close $cp;

	filter_BLAST();
}

sub filter_BLAST {
	my @linePart;
	my %unique;
	my $totalLines;
	my $currentLine;

	print "Filtering BLAST output.\n";
	# print "Opening files.\n";

	if ($idba_splited) {
		for (my $i = 1; $i <= $splitCount; $i++) {
			open(my $fh, "<", "$blast_output_dir/round$roundOfUnaligned/part$i/blastOUT") or die "Cannot open BLAST output.\n";
			open(my $parsed, ">", "$blast_output_dir/round$roundOfUnaligned/part$i/blastOUT_parsed.txt") or die "Cannot create BLAST filter output file.\n";

			while (<$fh>) {
				$totalLines++;
			}

			close $fh;

			open(my $fh2, "<", "$blast_output_dir/round$roundOfUnaligned/part$i/blastOUT") or die "Cannot open BLAST output.\n";

			$totalLines = int($totalLines / 2);
			$currentLine = 1;

			while (<$fh2>) {
				chomp;

				@linePart = split("\t");

				if (! $unique{$linePart[1]}) {
					print $parsed $linePart[1] . "\n";
					$unique{$linePart[1]}++;
				}

				if ($currentLine == $totalLines) {
					last;
				}

				$currentLine++;
			}

			close $fh2;
			close $parsed;

			system("$blast_location/blastdbcmd -db nt -entry_batch $blast_output_dir/round$roundOfUnaligned/part$i/blastOUT_parsed.txt -out $blast_output_dir/round$roundOfUnaligned/part$i/out.fna -outfmt %f");
		}

	} else {
		open (my $fh, "<", "$blast_output_dir/round$roundOfUnaligned/blastOUT") or die "Cannot open BLAST output.\n";
		open (my $parsed, ">", "$blast_output_dir/round$roundOfUnaligned/blastOUT_parsed.txt") or die "Cannot create BLAST filter output file.\n";

		# print "Counting lines.\n";

		while (<$fh>) {
			$totalLines++;
		}

		close $fh;

		# print "Opening fh2.\n";

		open (my $fh2, "<", "$blast_output_dir/round$roundOfUnaligned/blastOUT") or die "Cannot open BLAST output.\n";

		$totalLines = int($totalLines / 2);
		$currentLine = 1;

		while (<$fh2>) {
			# print "Loop 1.\n";
			chomp;

			@linePart = split("\t");

			if (! $unique{$linePart[1]}) {
				print $parsed $linePart[1] . "\n";
				$unique{$linePart[1]}++;
			}

			# print "Loop 2.\n";

			if ($currentLine == $totalLines) {
				last;
			}

			# print "Loop 3.\n";

			$currentLine++;
		}

		close $fh2;
		close $parsed;

		print "Retrieving fna sequence from BLAST database. ";

		system("$blast_location/blastdbcmd -db nt -entry_batch $blast_output_dir/round$roundOfUnaligned/blastOUT_parsed.txt -out $blast_output_dir/round$roundOfUnaligned/out.fna -outfmt %f");
		print "Done.\n";
	}

	parseBLASTfna();

	print "Done.\n";
}

sub parseBLASTfna {
	my $fileName;
	my $output;

	print "Parsing BLAST fna sequences. ";

	open(my $input, "<", "$blast_output_dir/round$roundOfUnaligned/out.fna");
	mkdir("$blast_output_dir/round$roundOfUnaligned/fnaFiles");

	for my $line (<$input>) {
		if (substr($line, 0, 1) eq ">") {
			$fileName = (split(/\./, (split(/\|/, $line))[3]))[0] . ".fna";

			if (defined $output) {
				if (fileno $output) {
					close $output;
				}
			}

			open($output, ">", "$blast_output_dir/round$roundOfUnaligned/fnaFiles/$fileName");
		}

		print $output $line;
	}

	print "Done.\n";
}

sub separateSAM {
	print "Separating SAM file into individual species' SAM file. ";

	unless (-d "$bowtie_output_dir/round$roundOfUnaligned/samSplit") {
		mkdir "$bowtie_output_dir/round$roundOfUnaligned/samSplit";
	}

	my @entry;

	open(SAMhandle, "<", $SAMfile);

	while (<SAMhandle>) {
		chomp;
		my @entry = split("\t");

		if (substr($entry[0], 0, 1) eq "@") {
			next;
		}

		if ($entry[2] eq "*") {
			next;
		}	

		my $fileName = (split(/\|/, $entry[2]))[3];
		# print "fileName 1: $fileName\n";
		$fileName = (split(/\./, $fileName))[0];
		# print "fileName 2: $fileName\n";

		open(my $output, ">>", "$bowtie_output_dir/round$roundOfUnaligned/samSplit/$fileName.sam");

		for (my $i = 0; $i < 12; $i++) {
			# print $output join("\t", @entry) . "\n";
			print $output $entry[$i];
			print $output "\t";
		}
		print $output "\n";

		close $output;
	}

	print "Done.\n";
}

sub processFna {
	# open(my $input, "<", "$blast_output_dir/round$roundOfUnaligned/fnaFiles");

	# while (<$input>) {
	# 	chomp;
	# 	my @parts = split("/");
	# 	my $ID = (split(/\./, $parts[-1]))[0];
	# 	# print "ID is $ID\nPausing.\n";
	# 	# <STDIN>;

	# 	$fnaHash{$ID} = join("/", @parts);
	# }
	if ($roundOfUnaligned == 1) {
		open(my $input, "<", "$kraken_output_dir/listOfFnaFiles");

		while (<$input>) {
			chomp;
			my @parts = split("/");
			my $ID = (split(/\./, $parts[-1]))[0];
			# print "ID is $ID\nPausing.\n";
			# <STDIN>;

			$fnaHash{$ID} = join("/", @parts);
		}
	} else {
		# open(my $input, "<", "$blast_output_dir/round" . ($roundOfUnaligned - 1) . "/fnaFiles");
		chdir "$blast_output_dir/round" . ($roundOfUnaligned - 1) . "/fnaFiles";

		my @files = <*.fna>;

		for my $file (@files) {
			chomp $file;
			my $ID = (split(/\./, $file))[0];

			# print "ID is $ID.";
			# print "fnaHash{ID} is $blast_output_dir/round$roundOfUnaligned/fnaFiles/$file.\nPausing.\n";
			# <STDIN>;

			$fnaHash{$ID} = "$blast_output_dir/round" . ($roundOfUnaligned - 1) . "/fnaFiles/$file";
		}
	}

	# print Dumper(\%fnaHash);

	chdir $home_dir;
}

sub runCompression {
	print "Running CRAM.\n";

	unless (-d $cram_dir) {
		mkdir $cram_dir or die "Output directory cannot be created.\n";
	}

	processFna();

	my @samFiles = <$bowtie_output_dir/round$roundOfUnaligned/samSplit/*.sam>;

	if (! @samFiles) {
		die "diamond brackets don't work.\n";
	}

	chdir $cram_dir;

	for my $file (@samFiles) {
		my $ID = (split(/\./, (split("/", $file))[-1]))[0];

		# print "ID is $ID.\nfnaHash{ID} is $fnaHash{$ID}.\ncwd is " . getcwd . "\nfile is $file\n";
		# print " Pausing.\n";
		# <STDIN>;

		# copy($fnaHash{$ID}, getcwd);
		system("cp $fnaHash{$ID} $cram_dir");
		# rename("$ID.fna", "$ID.fa");
		system("mv $cram_dir/$ID.fna $cram_dir/$ID.fa");

		system("samtools view -bT $cram_dir/$ID.fa $bowtie_output_dir/round$roundOfUnaligned/samSplit/$ID.sam 2> /dev/null | samtools sort - $cram_dir/Sorted$ID");
		# system("java -jar $cram cram --reference-fasta-file $cram_dir/$ID.fa --input-bam-file $cram_dir/Sorted$ID.bam --output-cram-file $cram_dir/$ID.cram");

		if ($huffman) {
			 system("java -jar $huffmanPath cram --capture-all-tags -Q -n -I $cram_dir/Sorted$ID.bam -R $cram_dir/$ID.fa -O $cram_dir/$ID.cram");
			# system("java -jar $huffmanPath cram -n -I $cram_dir/Sorted$ID.bam -R $cram_dir/$ID.fa -O $cram_dir/$ID.cram");
			system("touch .huffman");
		} elsif ($golomb) {
			system("java -jar $golombPath cram --capture-all-tags -Q -n -I $cram_dir/Sorted$ID.bam -R $cram_dir/$ID.fa -O $cram_dir/$ID.cram");
			system("touch .golomb");
		} elsif ($exGolomb) {
			system("java -jar $exGolombPath cram --capture-all-tags -Q -n -I $cram_dir/Sorted$ID.bam -R $cram_dir/$ID.fa -O $cram_dir/$ID.cram > $ID.txt");
			system("touch .exGolomb");
		}
	}

	print "Done.\n";

	$semaphore->down();

	print "Running MFCompress.\n";

	if ($roundOfUnaligned == $maxRounds) {
		system("cp $SAMfilter_dir/unaligned$roundOfUnaligned $cram_dir");
		system("$mfcPath $cram_dir/unaligned$roundOfUnaligned");
		system("rm $cram_dir/unaligned$roundOfUnaligned");
	}

	$semaphore->up();

	print "Done.\n";
}

sub doUnaligned {
	$semaphore->down();
	read_start_position();
	$semaphore->up();

	# if ($roundOfUnaligned != 1) {
	# 	# DEBUG: disable this
	# 	read_map_fromBlast();
	# } else {
	# 	# DEBUG: disable this
	# 	read_map();
	# }

	if ($roundOfUnaligned < $maxRounds) {
		run_IDBA();
		run_BLAST();
	}
}

sub processOutput {
	for (my $i = 1; $i <= $maxRounds; $i++) {
		mkdir("$output_dir_temp/Round$i") or die "Output directory cannot be created.\n";

		if ($i == 1) {
			system("cp $kraken_output_dir/listOfFnaFiles $output_dir_temp/Round$i");
		} else {
			system("cp $blast_output_dir/round" . ($i - 1) . "/blastOUT_parsed.txt $output_dir_temp/Round$i");
		}

		system("cp $output_dir_temp/cramOutput/round$i/*.cram $output_dir_temp/Round$i");

		if ($i == $maxRounds) {
			system("cp $output_dir_temp/cramOutput/round$i/*.mfc $output_dir_temp/Round$i");
		}

		if ($huffman) {
			system("cp $output_dir_temp/cramOutput/round$i/.huffman $output_dir_temp/Round$i");
		} elsif ($golomb) {
			system("cp $output_dir_temp/cramOutput/round$i/.golomb $output_dir_temp/Round$i");
		} elsif ($exGolomb) {
			system("cp $output_dir_temp/cramOutput/round$i/.exGolomb $output_dir_temp/Round$i");
		}
	}

    # TODO
    # do preservation of bowtie alignment rate, IDBA contigs, and bowtie samfiles
    # for alignment rate, maybe try redirecting using &> or 2>
    # http://stackoverflow.com/questions/418896/how-to-redirect-output-to-a-file-and-stdout
}

run_kraken();
filter_kraken();

while ($roundOfUnaligned <= $maxRounds) { 
	print "Round $roundOfUnaligned\n";
	run_bowtie();
	# DEBUG: disable this
	# read_start_position();

	# if ($roundOfUnaligned != 1) {
	# 	# DEBUG: disable this
	# 	read_map_fromBlast();
	# } else {
	# 	# DEBUG: disable this
	# 	read_map();
	# }

	# $R = Statistics::R->new();
	# $R -> start();

	# runR();
	# exit 0;
	# DEBUG: disable this
	# run_IDBA();
	# DEBUG: disable this
	# run_BLAST();
	# runCompression();
	# exit 0;
	threads->create(\&doUnaligned);
	threads->create(\&runCompression);

	foreach my $thread (threads->list()) {
		$thread -> join();
	}

	$roundOfUnaligned++;
}

processOutput();

exit 0;

sub runR {
	my $R;

	print "Running R. ";

	chdir($SAMfilter_dir);
	my @startposFiles = read_dir($SAMfilter_dir, chomp => 1);
	open(my $distFile, ">>", "startpos.distribution");
	# $R = Statistics::R->new();
	# $R -> start();
 
	# if (! $R->is_started()) {
		# print "R did not start. Trying again.\n";
		# $R -> start();
	# }

	for my $startposFile (@startposFiles) {
		if (-d $startposFile) {
			next;
		}

		my $extension = (split(/\./, $startposFile))[1];

		if (! defined $extension) {
			next;
		}

		if ($extension ne "startpos" || (split(/\./, $startposFile))[2]) {
			next;
		}

		if ((-s $startposFile) < 3000) {
			next;
		}

		$R = Statistics::R->new();

		# print "File is: $startposFile\n";
		
		# open(my $distFile, ">", "$startposFile.distribution");

		$R -> set('INPUT', $startposFile);
		# $R -> run($cmds);
		my $out = $R -> run_from_file("/shared3/readStartPosDistribution.R");
		
		# print "\n" . $out . "\n";
		# print "PAUSED.\n";
		# <STDIN>;

		my $bestDist = $R -> get('bestDist');
		my $errorCode = $R -> get('errorCode');
		my $dataTest_ref = $R -> get('dataTest');
		# my $ks1 = $R -> get('ks1');
		# my $ks2 = $R -> get('ks2');
		# my $ks3 = $R -> get('ks3');
		# my $fnbinom = $R -> get('fnbinom');
		# my $fgeom = $R -> get('fgeom');
		# my $fpois = $R -> get('fpois');
		# my $dataTest = $R -> get('dataTest');

		# print "File: $startposFile\n";

		# if (ref($bestDist) eq "ARRAY") {
		# 	print "bestDist is: ";
		# 	print join("," , @{ $bestDist });
		# 	print "\n";
		# } else {
		# 	print "bestDist is $bestDist.\n";
		# }

		# print "R errorCode is $errorCode.\n";
		# print "dataTest is: @$dataTest_ref\n";

		# if (ref($ks1) eq "ARRAY") {
		# 	print "ks1 is: ";
		# 	print join("," , @{ $ks1 });
		# 	print "\n";
		# } else {
		# 	print "KS1 = " . $R -> get('ks1$statistic[[1]]') . "\n";
		# }

		# if (ref($ks2) eq "ARRAY") {
		# 	print "ks2 is: ";
		# 	print join("," , @{ $ks2 });
		# 	print "\n";
		# } else {
		# 	print "KS2 = " . $R -> get('ks2$statistic[[1]]') . "\n";
		# }

		# if (ref($ks3) eq "ARRAY") {
		# 	print "ks3 is: ";
		# 	print join("," , @{ $ks3 });
		# 	print "\n";
		# } else {
		# 	print "KS3 = " . $R -> get('ks3$statistic[[1]]') . "\n";
		# }

		# if (ref($fnbinom) eq "ARRAY") {
		# 	print "fnbinom is: ";
		# 	print join("," , @{ $fnbinom });
		# 	print "\n";
		# } else {
		# 	print "fnbinom = " . $R -> get('fnbinom[1]') . "\n";
		# }

		# if (ref($fgeom) eq "ARRAY") {
		# 	print "fgeom is: ";
		# 	print join("," , @{ $fgeom });
		# 	print "\n";
		# } else {
		# 	print "fgeom = " . $R -> get('fgeom[1]') . "\n";
		# }

		# if (ref($fpois) eq "ARRAY") {
		# 	print "fpois is: ";
		# 	print join("," , @{ $fpois });
		# 	print "\n";
		# } else {
		# 	print "fpois = " . $R -> get('fpois[1]') . "\n";
		# }

		# if (ref($dataTest) eq "ARRAY") {
		# 	print "dataTest is: ";
		# 	print join("," , @{ $dataTest });
		# 	print "\n";
		# } else {
		# 	print "dataTest = " . $R -> get('dataTest') . "\n";
		# }

		# <STDIN>;

		print $distFile ">$startposFile\n";

		if ($bestDist == 1) {
			# print $R -> get('method') . "distribution.\n";
			# print "r is " . $R -> get('r') . "\n";
			# print "p is " . $R -> get('p') . "\n\n";
			$startposHash{$startposFile} = join("\t", $R -> get('method'), $R -> get('r'), $R -> get('p'), @$dataTest_ref);
			print $distFile $startposHash{$startposFile};
			# close $distFile;
		} elsif ($bestDist == 2) {
			# print $R -> get('method') . "distribution.\n";
			# print "p is " . $R -> get('p') . "\n\n";
			$startposHash{$startposFile} = join("\t", $R -> get('method'), $R -> get('p'), @$dataTest_ref);
			print $distFile $startposHash{$startposFile};
			# close $distFile;
		} elsif ($bestDist == 3) {
			# print $R -> get('method') . "distribution.\n";
			# print "lambda is " . $R -> get('lambda') . "\n\n";
			$startposHash{$startposFile} = join("\t", $R -> get('method'), $R -> get('lambda'), @$dataTest_ref);
			print $distFile $startposHash{$startposFile};
			# close $distFile;
		} elsif ($bestDist == -1) {
			# print "File frequency is " . $R -> get('errorLength') . "\n";
			# $startposHash{$startposFile} = join("\t", $R -> get('method'), $R -> get('lambda'));
			# print $distFile $startposHash{$startposFile};
			# close $distFile;
		} else {
			print "*** R did not return any distribution on file $startposFile. ***\n";
		}

		if ($errorCode == 1) {
			# print "Negative binomial test failed.\n";
			print $distFile "\nNegative binomial test failed.\n";
		} else {
			# print "Negative binomial test passed.\n";
			print $distFile "\nNegative binomial test passed.\n"
		}

		# close $distFile;
		$R -> stop();
	}

	# $R->stop();
	close $distFile;

	print "Done.\n";
}

sub read_map {
	# DEBUG: disable this
	unless (-d $readMap_dir) {
		mkdir($readMap_dir) or die $!;
	}

	if (-e "$readMap_dir/.cpRM") {
		print "Read map checkpoint file already exists. Skipping.\n";
		return;
	}

	print "Read mapping. ";

	# Replacement for read_mapper_compression.sh

	# TODO this part is easily multi-threadable by converting the foreach loop into a sub and divide SAMfilteredFiles by threads

	@SAMfilteredFiles = read_dir($SAMfilter_dir, chomp => 1);
	@genus = read_dir($genus_db_location, chomp => 1, prefix => 1);

	foreach my $currentFile (@SAMfilteredFiles) {
		my @splitName = split(/\./, $currentFile);

		if (! $splitName[1]) {
			next;
		}

		if ($splitName[1] ne "alignments") {
			#print("Skipping $currentFile\n");
			next;
		}

		$fileID = $splitName[0];

		# print("Finding $currentFile. fileID is $fileID\n");
		# <STDIN>;

		eval {
			find(\&wanted, @genus);
		};
	}

	print "Done.\n";

	open(my $cpRM, ">", "$readMap_dir/.cpRM") or die "Read map checkpoint file cannot be created.\n";
	close $cpRM;
}

# this is read_mapper_compression.sh
sub wanted {
	if (-d) {
		return;
	}

	# print("current file is: $File::Find::name\n");
	my $ID = (split(/\./))[0];
	# print("current ID is: $ID\n");

	if ($ID eq $fileID) {
		#change this to the SAMfilter directory
		my $genusName = (split("_", (split("/", $File::Find::dir))[-1]))[0];
		# print("Genus name is: $genusName\n");
		# print("\nMatched $_ to $fileID\n");

		open(my $fh, ">>", "$readMap_dir/$genusName") or die $!;
		# print("Making/printing to $SAMfilter_dir/$genusName\n");
		# print "$ID.alignments\n";
		# <STDIN>;

		print $fh "$ID.alignments\n";
		close $fh;
		# push(@madeFiles, $ID);

		open(my $fh2, ">>", "$readMap_dir/" . $genusName . "_combined.fna");
		open(INPUT, "$SAMfilter_dir/$ID.alignments");

		while (<INPUT>) {
			chomp;
			print $fh2 "$_\n";
		}

		close $fh2;
		die;
	}
}

sub read_map_fromBlast {
	getHeaders() or die $!;
	# print "Size: @headers\n";
	# exit 0;
	$readMap_dir = "$SAMfilter_dir/readMap";
	# DEBUG: disable this
	unless (-d $readMap_dir) {
		mkdir($readMap_dir);
	}

	print "Read mapping. ";

	@SAMfilteredFiles = read_dir($SAMfilter_dir, chomp => 1);

	foreach my $currentFile (@SAMfilteredFiles) {
		my @splitName = split(/\./, $currentFile);
		# print("Currentfile: $currentFile\n");

		if (! $splitName[1]) {
			# print("Skipping $currentFile\n");
			next;
		}

		if ($splitName[1] ne "alignments") {
			# print("Skipping $currentFile\n");
			next;
		}

		# <STDIN>;

		$fileID = $splitName[0];

		# print("Finding $currentFile. fileID is $fileID\n");
		# <STDIN>;

		for my $header (@headers) {
			my @headerParts = split(" ", $header);
			my $ID = (split(/\|/, $headerParts[0]))[3];
			$ID = (split(/\./, $ID))[0];

			if ($ID eq $fileID) {
				# print "$fileID matches $ID\n";
				# <STDIN>;

				open(my $fh, ">>", "$readMap_dir/$headerParts[1]") or die $!;

				print $fh "$ID.alignments\n";
				close $fh;
				# push(@madeFiles, $ID);

				open(my $fh2, ">>", "$readMap_dir/$headerParts[1]_combined.fna") or die $!;
				open(INPUT, "$SAMfilter_dir/$ID.alignments") or die $!;

				while (<INPUT>) {
					chomp;
					print $fh2 "$_\n";
				}

				close $fh2;

				last;
			}
		}
	}

	print "Done.\n";
}

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub do_decompress {
	printf "decompressing is still under construction :)\n";
	exit 0;

	for (my $i = 1; $i <= $maxRounds; $i++) {
		# unless (-d)
	}
}

sub decompressCRAM {
	for (my $i = 1; $i <= $maxRounds; $i++) {

	}
}