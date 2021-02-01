use strict;
use warnings;
use Cwd;

my $gb_file = "NC_045512.2.gb";
my $fasta_input = "210128_nucleotides.fasta";

my $outdir = "input";

my $pwd = getcwd();

my $seqs = ParseFasta($fasta_input);

mkdir $outdir;
chdir $outdir;

foreach my $id (keys %$seqs) {
	$id =~ m/^(\S+)/i;
	my $outfile = $1;
	$outfile .= ".fas";
	
	open (my $out_fh, ">", $outfile) or die $!;
	print $out_fh ">$id\n", $seqs -> {$id}, "\n";
	close $out_fh;
}

chdir $pwd;

sub ParseFasta {
	my $file = shift;
	
	my $sequences;
	
	open (my $in_fh, "<", $file) or die $!;
	my $name;
	while (my $line = <$in_fh>) {
		chomp $line;
		$line =~ s/\R+//g;		
		if ($line =~ m/^>(.+)/i) {
			$name = $1;
		}
		else {
			$sequences -> {$name} .= $line;
		}
	}
	close $in_fh;
	return $sequences;
}