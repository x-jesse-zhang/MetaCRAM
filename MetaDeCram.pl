#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $inputDir;
my $blast_location = "/opt/ncbi-blast-2.2.29+/bin";
my %fnaHash;
my $huffmanPath = getcwd . "/HuffmanCRAM.jar";
my $golombPath = getcwd . "/GolombCRAM.jar";
my $exGolombPath = getcwd . "/ExtendedGolombCRAM.jar";
my $mfdPath = getcwd . "/MFCompressD";

my $maxRounds = 2;
my $fileCounter = 0;

my $DEBUGLOGGING = 1;

GetOptions(
	"input=s" => \$inputDir
);

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
} 

chdir $inputDir;
open(my $output1, ">", "OUT");
open(my $output2, ">", "OUT2");

for (my $i = 1; $i <= $maxRounds; $i++) {
	print "ROUND $i\n";
	# print STDERR "ROUND $i\n";

	# if ($i == 2) {
	# 	print STDERR "PAUSING\n";
	# 	<STDIN>;
	# 	print STDERR "Unpaused\n";
	# }
	# <STDIN>;
	# print "Unpaused\n";

	my $currentDir = $inputDir . "/Round$i";
	chdir $currentDir;
	# print STDERR "currentDir: $currentDir\n";

	# fasta file processing
	if ($i == 1) {
		print "Acquiring fasta files from Kraken\n";

		if (! -e "listOfFnaFiles") {
			print "Cannot acquire reference fasta files for CRAM decompression round$i at $currentDir.\n";
			exit(1);
		}

		if (! -e ".huffman" && ! -e ".golomb" && ! -e ".exGolomb") {
			print "Cannot find compression method for CRAM decompression round$i at $currentDir.\n";
			exit(2);
		}

		open(my $input, "<", "listOfFnaFiles");

		while(<$input>) {
			$fileCounter++;
			chomp;
			my @parts = split("/");
			my $ID = (split(/\./, $parts[-1]))[0];

			$fnaHash{$ID} = join("/", @parts);
		}

		close $input;

		# print "Acquired $fileCounter fasta files\n";
	} else {
		print "Acquiring fasta files from blastdbcmd\n";

		if (! -e "blastOUT_parsed.txt") {
			print "Cannot acquire reference fasta files for CRAM decompression round$i at $currentDir.\n";
			exit(1);
		}

		if (! -e ".huffman" && ! -e ".golomb" && ! -e ".exGolomb") {
			print "Cannot find compression method for CRAM decompression round$i at $currentDir.\n";
			exit(2);
		}

		chdir "/shared3";
		system("$blast_location/blastdbcmd -db nt -entry_batch $currentDir/blastOUT_parsed.txt -out $currentDir/out.fna -outfmt %f");
		chdir $currentDir;
		# print "Pausing after blastdbcmd\n";
		# <STDIN>;
		# print "Unpaused\n";

		my $fileName;
		my $fastaFile;

		open(my $input, "<", "out.fna");
		# mkdir("fasta");

		for my $line (<$input>) {
			if (substr($line, 0, 1) eq ">") {
				$fileName = (split(/\./, (split(/\|/, $line))[3]))[0] . ".fa";
				# print "Pausing. filename: " . $fileName . "\n";
				# <STDIN>;
				# print "Unpause";

				if (defined $fastaFile) {
					if (fileno $fastaFile) {
						close $fastaFile;
					}
				}

				open($fastaFile, ">", "$fileName");
			}

			print $fastaFile $line;
		}

		close $fastaFile;

		my @files = <*.fa>;

		# print "Acquired " . scalar @files . " fasta files.\n";

		for my $file (@files) {
			chomp $file;
			my $ID = (split(/\./, $file))[0];

			# print "ID is $ID.";
			# print "fnaHash{ID} is $file.\nPausing.\n";
			# <STDIN>;
			# print "Unpaused.\n";

			$fnaHash{$ID} = $file;
		}
	}

	my @cramFiles = <*.cram>;

	print "CRAM decompressing. " . scalar @cramFiles . " cram files present.\n";
	# print STDERR "CRAM decompressing. " . scalar @cramFiles . " cram files present.\n";

	for my $file (@cramFiles) {
		my $ID = (split(/\./, (split("/", $file))[-1]))[0];

		if ($i == 1) {
			system("cp $fnaHash{$ID} $currentDir");

			if (-e "$ID.fna") {
				system("mv $ID.fna $ID.fa");
			}
		}

		# print "ID: " . $ID . "\n";
		# <STDIN>;

		system("samtools faidx $ID.fa");

		if (-e ".huffman") {
			system("java -jar $huffmanPath bam --input-cram-file $ID.cram --reference-fasta-file $ID.fa --print-sam-header > $ID.sam");
		} elsif (-e ".golomb") {
			system("java -jar $golombPath bam --input-cram-file $ID.cram --reference-fasta-file $ID.fa --print-sam-header > $ID.sam");
		} elsif (-e ".exGolomb") {
			system("java -jar $exGolombPath bam --input-cram-file $ID.cram --reference-fasta-file $ID.fa --print-sam-header > $ID.sam");
		} else {
			die "error.\n";
		}

		open(my $sam, "<", "$ID.sam");
		$fileCounter = 0;

		print "On $ID.sam. Round$i. \n";

		while (my $line = <$sam>) {
			chomp $line;
			my @lineParts = split("\t", $line);

			if (substr($lineParts[0], 0, 1) eq "@") {
				next;
			}

			# if (!$output1->opened()) {
			# 	print STDERR "output1 is closed.\n";
			# 	die;
			# } elsif (!$output2->opened()) {
			# 	print STDERR "output2 is closed.\n";
			# 	die;
			# }

			select((select($output1), $| = 1)[0]);
			select((select($output2), $| = 1)[0]);

			# if current mate is mate 2 by ID or by matching
			# if (($lineParts[1] & 128) == 128 || index($lineParts[0], "/2") != -1) {
			if (($lineParts[1] & 128) == 128) {
				# select(STDOUT);
				if ($DEBUGLOGGING) {
					print "Printing to OUT2 for Round$i from $ID.sam: $lineParts[0]\n";
				}

				# select($output2);
				print $output2 ">$lineParts[0]\n" or die "printing failed\n";
				# print $output2 "$lineParts[9]\n" or die "printing failed\n";
				# 16 is the tag for aligning to reverse strand
				if (($lineParts[1] & 16) == 16) {
					print $output2 reverse_complement($lineParts[9]) . "\n" or die "printing failed\n";
				} else {
					print $output2 "$lineParts[9]\n" or die "printing failed\n";
				}
			} else {	# otherwise current mate is mate 1
				# select(STDOUT);
				if ($DEBUGLOGGING) {
					print "Printing to OUT for Round$i from $ID.sam: $lineParts[0]\n";
				}

				# select($output1);
				if (index($lineParts[0], "/0") != -1) {
					print $output1 ">" . substr($lineParts[0], 0, index($lineParts[0], "/0")) . "\n";
				} else {
					print $output1 ">$lineParts[0]\n";
				}
				
				if (($lineParts[1] & 16) == 16) {
					print $output1 reverse_complement($lineParts[9]) . "\n" or die "printing failed\n";
				} else {
					print $output1 "$lineParts[9]\n" or die "printing failed\n";
				}
			}

			$fileCounter++;	
			# select(STDOUT);		
		}

		close $sam;

		# print "$fileCounter lines printed for $ID.sam\n";
	}

	# unaligned file processing
	if ($i == 2) {
		print "Processing unaligned file.\n";
		# print STDERR "Processing unaligned file.\n";

		my @mfcFiles = <*.mfc>;

		for my $file (@mfcFiles) {
			system("$mfdPath $file");
			my $originalFile = $file . ".d";
			# print "file: $file originalFile: $originalFile. Pausing.\n";
			# <STDIN>;
			# print "Unpaused\n";

			open(my $unaligned, "<", $originalFile);

			while (my $line = <$unaligned>) {
				chomp $line;
				my @lineParts = split("\t", $line);

				# if (index($lineParts[0], "/1") != 1) {
				# 	# select($output2);
				# 	print $output2 "$line\n";
				# 	$line = <$unaligned>;
				# 	print $output2 "$line\n";
				# } else {
				# 	# select($output1);
				# 	print $output1 "$line\n";
				# 	$line = <$unaligned>;
				# 	print $output1 "$line\n";
				# }

				if (index($lineParts[0], "/0") == -1) {
					# select($output2);
					if ($DEBUGLOGGING) {
						print "Printing to OUT2 for Round$i for unaligned: $line\n";
					}
					print $output2 "$line\n";
					$line = <$unaligned>;
					chomp $line;
					print $output2 "$line\n";
				} else {
					# select($output1);
					# print $output1 "$line\n";
					if ($DEBUGLOGGING) {
						print "Printing to OUT2 for Round$i for unaligned: $line\n";
					}
					print $output1 substr($line, 0, index($line, "/0")) . "\n";
					$line = <$unaligned>;
					chomp $line;
					print $output1 "$line\n";
				}

				# select(STDOUT);
			}

			close $unaligned;	
		}	
	}
}

close $output1;
close $output2;

if ($DEBUGLOGGING) {
	system("perl /shared3/OUTsorter.pl $inputDir/OUT $inputDir/sortedOUT");
	system("perl /shared3/OUTsorter.pl $inputDir/OUT2 $inputDir/sortedOUT2");
}
# ideally now we have SAM files and the unaligned files. now to create parser...
