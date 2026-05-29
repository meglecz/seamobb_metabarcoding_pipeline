use warnings;
use strict;
use Data::Dumper;


# AIM:
# Pool different ASV tables into one
# Cluster ASV using eiter of the two methods 
#		-vsearch --cluster_fast using fixed identity threshold
#		-SWARM 
# Add a custerid and custersize (number of ASV in the cluster) column to the pooled ASV table

# NOTES:
# If same sample name in different runs, keep them separatelly
# No taxassign columns in the input files
# Last column is the sequence column in input files
# Works only for ASV tables of the same marker
# Include a column with the list of runs, where the sample is present

# INPUT:
# ASV tables
#  Output of vtam filter
#  Columns: run	marker	variant	sequence_length	read_count	(list_of_samples)
#			clusterid	clustersize	chimera_borderline	sequence


# Third party prg:
# 	vsearch if $cluster_method = 'cluster_fast'
# 	swarm if $cluster_method = 'swarm'

my %files = (
	'Arm01' => '/home/meglecz/seamobb_metabarcoding_pipeline/results/Arm01/5_filter_optimized/asvtable.tsv',
	'Arm02' => '/home/meglecz/seamobb_metabarcoding_pipeline/results/Arm02/5_filter_optimized/asvtable.tsv',
	'Arm03' => '/home/meglecz/seamobb_metabarcoding_pipeline/results/Arm03/5_filter_optimized/asvtable.tsv',
	'Arm04' => '/home/meglecz/seamobb_metabarcoding_pipeline/results/Arm04/5_filter_optimized/asvtable.tsv',
	'TEST1' => '/home/meglecz/seamobb_metabarcoding_pipeline/results/TEST1/5_filter_optimized/asvtable.tsv',
);
my $outfile = "/home/meglecz/seamobb_metabarcoding_pipeline/results/pooled_asv_tables.tsv";
my $tmp_dir = "/home/meglecz/seamobb_metabarcoding_pipeline/results/tmp/";
my $marker = 'IIICBR';
my $cluster_method = 'swarm'; #swarm/cluster_fast
my $cluster_limit = 7; # [0-1]if cluster_fast, positive integer if swarm

my $first_sample = 5;
my $after_last_sample = '^clusterid';
my $sid = 2;


##### read data from input asvtables
my %hash; # hash{var_seq}{run-sample} = readcount
my %seq_id; #$seq_id{seq} = id;
my %run_samples; #$run_samples{run-sample} = ''
my %varid_run; #%varid_run{varid}{run} = number of sample where the ASV is present
my %sample_count;  #%varid_run{varid}{run-sample} = ''
my %rc; # $rc{varid} = total number of reads;

foreach my $run (sort keys %files)
{
	my $file = $files{$run};
	my %title; # title{sample} = run-sample
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>;
	$title =~ s/\s*$//;
	my @title = split("\t", $title);
	my $last_sample_ind = get_column_index($title, $after_last_sample, "\t") -1;

	for(my $i = $first_sample; $i<= $last_sample_ind; ++$i )
	{
		my $run_sample = $run.'-'.$title[$i];
		$run_samples{$run_sample} = '';
		$title{$title[$i]} = $run_sample;
	}
	
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		$seq_id{$line[-1]} = $line[$sid];
		for(my $i = $first_sample; $i <= $last_sample_ind; ++$i)
		{
			if($line[$i])
			{
				$hash{$line[-1]}{$title{$title[$i]}} = $line[$i]; # $hash{seq}{run-sample} = readcount
				$varid_run{$line[$sid]}{$run} = '';
				$sample_count{$line[$sid]}{$title{$title[$i]}} = '';
				$rc{$line[$sid]} += $line[$i];
			}
		}
	}
	close IN;
}
#print Dumper(\%sample_count);


##### clustering
my %cluster_info; # $cluster_info{varid} = (clusterid, clustersize)
if($cluster_method eq 'cluster_fast')
{
	%cluster_info = add_cluster($tmp_dir, $cluster_limit, \%seq_id);
}
elsif($cluster_method eq 'swarm')
{
	my $c = 'mkdir '.$tmp_dir;
	system $c;
	
	my $fas = $tmp_dir.'swarm_input.fas';
	open(FAS, '>', $fas) or die "Cannot open $fas\n";
	foreach my $sequence (keys %seq_id)
	{
		my $sid = $seq_id{$sequence};
		print FAS ">$sid", '_', $rc{$sid} , "\n$sequence\n";
	}
	close FAS;
	
	my $swarm_o = $tmp_dir.'swarm_output_'.$cluster_limit.'.txt';
	my $cmd = 'swarm -t 8 -d '.$cluster_limit.' -o '.$swarm_o;
	if($cluster_limit == 1)
	{
		$cmd .= ' -f';
	}
	$cmd .= ' '.$fas;
	system $cmd;
	read_swarm_output(\%cluster_info, $swarm_o);
#	print Dumper(\%cluster_info );
}
else 
{
	print "$cluster_method is not implemented\n";
	exit;
}


#### print pooled asv table inclusding clustering info
open(OUT, '>', $outfile) or die "Cannot open $outfile\n";
my @run_samples = sort keys %run_samples;
print OUT "marker	variant	sequence_length	read_count	sample_count	runs	", join("\t", @run_samples), "\tclusterid_",$cluster_limit,"\tclustersize	sequence\n";


foreach my $seq (sort keys %hash)
{
	my $rc = 0;
	my @temp = ();
	foreach my $run_sample (@run_samples)
	{
		if(exists $hash{$seq}{$run_sample})
		{
			push(@temp, $hash{$seq}{$run_sample});
			$rc += $hash{$seq}{$run_sample};
		}
		else
		{
			push(@temp, 0);
		}
	}
	my $varid = $seq_id{$seq};
	print OUT "$marker	$varid	", length $seq, "\t", $rc, "\t";
	print OUT scalar keys %{$sample_count{$varid}}, "\t", join(';', keys %{$varid_run{$varid}}), "\t";
	print OUT join("\t", @temp), "\t", join("\t", @{$cluster_info{$varid}}), "\t", "$seq\n";
}
close OUT;


exit;

######################################################
sub read_swarm_output
{
	my ($hash, $file) = @_;
my %hash; # $cluster_info{varid} = (clusterid, clustersize)
	open(IN, $file) or die "Cannot open $file\n";
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split(' ', $line);
		my $n = scalar @line;
		my $cent = shift @line; # get centroid
		$cent =~ s/_[0-9]+$//g; # delete read count
		@{$$hash{$cent}} = ($cent, $n);

		foreach my $var (@line)
		{
			$var =~ s/_[0-9]+$//g; 
			@{$$hash{$var}} = ($cent, $n);
		}
	}
	close IN;
}

######################################################

sub read_clusters
{
 my ($dir, $seq_clust, $centroids, $clusters) = @_;
 
#my %seq_clust; # $seq_clust{var_id} = cluster
#my %centroids; # $clust{clust_id}{centroid_id/centroidseq} = $centroid_id/$centroidseq
#my %clusters; # $clusters_varlist{clust_id}{varid} = '';

	my @files = get_file_list_from_folder($dir, '^cluster');
	
	foreach my $file (@files)
	{
		my $cl = $file;
		$file = $dir.$file;
		open(IN, $file) or die "Cannot open $file\n";
		my $centroid_id = <IN>; # get centroid ID
		$centroid_id =~ s/>//;
		$centroid_id =~ s/\s*$//;
		close IN;

		my %seq = read_fasta_to_hash_wo_space_gb($file); # read variants in the cluster
		
		$$centroids{$cl}{centroid_id} = $centroid_id;
		$$centroids{$cl}{centroidseq} = $seq{$centroid_id};

		foreach my $varid (keys %seq)  # fill $seq_clust{var_id} = cluster
		{
			$$seq_clust{$varid} = $cl;
			$$clusters{$cl}{$varid} = '';
		}
	}
}

#####################################################
sub read_fasta_to_hash_wo_space_gb
{
my ($filename) = @_;
my %seq = ();
	my $i = 0;
	open(IN, $filename) or die "cannot open $filename\n";
	$/ = ">";
	while (my $seq = <IN>)
	{
		$seq =~ s/>//;
		unless ($seq eq '')
		{
			$seq =~ s/.*\n//;
			my $code = $&;
#			$seq =~ s/$code//;
			$seq =~ s/\s//g;
			$code =~ s/\s.*//;
			$code =~ s/\s//;
			if (exists $seq{$code})
			{
				++$i;
#				print "$code\n";
			}
			$seq{$code} = $seq;
		}
	}
close IN;
	if ($i>0)
	{
		print "number of sequences with code already used by other sequence: $i\n";
	}
	$/ = "\n";

return %seq;
}

############################################

sub add_cluster
{
	my ($dir, $id, $seq) = @_;
	
# delete existing cluster folder to avoid mixing up data with previous runs
	my $temp = $dir.'tmp/';
	if (-e $temp) 
	{
		my $c = 'rm -r '.$temp;
		system $c;
	}
# make cluster folder
	my $c = 'mkdir '.$temp;
	system $c;

#Make a clustering with $id identity from all variants of both markers
	my $fas = $temp.'seq.fas';
	open(FAS, '>', $fas) or die "Cannot open $fas\n";
	foreach my $sequence (keys %$seq)
	{
		print FAS ">$$seq{$sequence}\n$sequence\n";
	}
	close FAS;
#	write_hash_to_fasta($fas, \%seq);
	my $cluster = 'vsearch --cluster_fast '.$fas.' --clusters '.$temp.'cluster --id '.$id;
	system $cluster;

# read cluster info to hashes
	my %seqid_clust; # $seq_clust{var_id} = cluster
	my %centroids; # $centroids{clust_id}{centroid_id/centroidseq} = $centroid_id/$centroidseq
	my %clusters_varlist; # $clusters_varlist{clust_id}{varid} = '';
	read_clusters($temp, \%seqid_clust, \%centroids, \%clusters_varlist);

	my %hash; #$hash{$seqid} = (clusterid, cluster_size)
	foreach my $seqid (keys %seqid_clust)
	{
		my $size = scalar keys  %{$clusters_varlist{$seqid_clust{$seqid}}};
		@{$hash{$seqid}} = ($seqid_clust{$seqid}, $size);
	}
	
	$c = 'rm -r '.$temp;
	system $c;

	
	return %hash;
}


##############################################
sub get_column_index
{
	my ($title, $motif, $sep) = @_;
	
	my @t = split($sep, $title);
	
	for(my $i = 0; $i < scalar @t; ++$i)
	{
		if($t[$i] =~ /$motif/)
		{
			return $i;
		}
	}

	print "$motif not found in the title line\n";
	exit;
}
##############################################
sub get_file_list_from_folder
{
 my ($folder, $file_motif) = @_;
 
  unless ( opendir(FOLDER, $folder) )
  {
      print "Cannot access to folder $folder\n";
      exit;
  }

my @filenames = grep ( !/^\.\.?$/, readdir(FOLDER) );

#print "@filenames\n";
closedir(FOLDER);
my @files = ();
foreach my $file (sort @filenames)
{
	if ($file =~ /$file_motif/)
	{
		push(@files, $file);
	}
}
@filenames = ();

#print "@files\n";
return @files;
}

