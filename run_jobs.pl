use strict;
use warnings;
use Cwd;
use File::Find;

my $gb_file = "NC_045512.2.gb";
my $search_path = "input";
my $array_file = "control_file.txt";
my $max_concurrent = 128;
my $script = "extract_seqs_array.sh";

open (my $out_fh, ">", $array_file) or die $!;

my $pwd = getcwd();
my $count = 0;

find (\&FindFasta, $search_path);

my $cmd_string = "qsub -tc $max_concurrent -t 1-$count:1 $script $array_file";
print "$cmd_string\n";
`$cmd_string`;

sub FindFasta {
    my $file = $_;
    my $full_path = $File::Find::name;
   
    if ($file =~ m/(.+)*\.fas/) {	    
	    my $name = $1;
	    my $output_file = "$pwd/output/$name.fas";
	    unless (-e $output_file) { 
	    	print $out_fh " $gb_file $full_path $name\n";
	    	$count++;
    	}
    }
}