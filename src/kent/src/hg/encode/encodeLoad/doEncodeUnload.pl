#!/usr/bin/env perl

# encodeUnload.pl - unload ENCODE data submission generated by the
#                       automated submission pipeline
# Reads load.ra for information about what to do

# Writes error or log information to STDOUT
# Returns 0 if unload succeeds.

# DO NOT EDIT the /cluster/bin/scripts copy of this file --
# edit the CVS'ed source at: ~/kent/src/hg/encode/encodeUnload/doEncodeUnload.pl
#
# $Id: doEncodeUnload.pl,v 1.8 2010/04/22 21:01:54 tdreszer Exp $

use warnings;
use strict;

use Getopt::Long;
use Cwd;
use File::Basename;

use FindBin qw($Bin);
use lib "$Bin";
use Encode;
use HgAutomate;
use HgDb;
use RAFile;


use vars qw/$opt_verbose $opt_configDir/;
my $PROG = basename $0;

sub usage
{
    print STDERR <<END;
usage: doEncodeUnload.pl pipeline-instance project_submission_dir

The pipeline instance variable is a switch that changes the behavior of doEncodeUnload.
The changes if the instance is:

standard
    allows use of hg19 and mm9 databases only

anything else
    allows use of the encodeTest database only

	project_submission_dir needs a full path
	OPTIONS:
		-verbose=i	Verbosity level
		-configDir=s	Config directory location
END
    exit(1);
}

sub genericUnload
{
    my ($assembly, $db, $tableName) = @_;
    $db->dropTableIfExist($tableName);
}

sub unloadWig
{
    my ($assembly, $db, $tableName) = @_;
    $db->dropTableIfExist($tableName);

    # remove symlink
    my $file = "/gbdb/$assembly/wib/$tableName.wib";
    if(-e $file) {
        HgAutomate::verbose(3, "removing wib '$file'\n");
        if(system("rm -f $file")) {
            die "unexpected error removing symlink $file";
        }
    }
}

sub unloadBam
{
    my ($assembly, $db, $tableName) = @_;
    $db->dropTableIfExist($tableName);

    # remove symlink
    my $file = "/gbdb/$assembly/bbi/$tableName.bam";
    if(-e $file) {
        HgAutomate::verbose(3, "removing bam '$file'\n");
        if(system("rm -f $file")) {
            die "unexpected error removing symlink $file";
        }
    }

    $file = "/gbdb/$assembly/bbi/$tableName.bam.bai";
    if(-e $file) {
        HgAutomate::verbose(3, "removing bai '$file'\n");
        if(system("rm -f $file")) {
            die "unexpected error removing symlink $file";
        }
    }
}

sub unloadBigWig
{
    my ($assembly, $db, $tableName) = @_;
    $db->dropTableIfExist($tableName);

    # remove symlink
    my $file = "/gbdb/$assembly/bbi/$tableName.bigWig";
    if(-e $file) {
        HgAutomate::verbose(3, "removing bigWig '$file'\n");
        if(system("rm -f $file")) {
            die "unexpected error removing symlink $file";
        }
    }
    # FIXME: Souldn't we remove files from downloads dir (and gbdb subdir) as well??
    #my $file = "/usr/local/apache/htdocs/goldenPath/$assembly/encodeDCC/encSydhTfbsStanf/gbdb/$tableName.bw";
    #if(-e $file) {
    #    HgAutomate::verbose(3, "removing wib '$file'\n");
    #    if(system("rm -f $file")) {
    #        die "unexpected error removing symlink $file";
    #    }
    #}
}

sub unloadBigBed
{
    my ($assembly, $db, $tableName) = @_;
    $db->dropTableIfExist($tableName);
    
    # remove symlink
    my $file = "/gbdb/$assembly/bbi/$tableName.bigBed";
    if(-e $file) {
        HgAutomate::verbose(3, "removing bigBed '$file'\n");
        if(system("rm -f $file")) { 
            die "unexpected error removing symlink $file";
        }
    }
    # FIXME: Souldn't we remove files from downloads dir (and gbdb subdir) as well??
    #my $file = "/usr/local/apache/htdocs/goldenPath/$assembly/encodeDCC/encSydhTfbsStanf/gbdb/$tableName.bw";
    #if(-e $file) {
    #    HgAutomate::verbose(3, "removing wib '$file'\n");
    #    if(system("rm -f $file")) {
    #        die "unexpected error removing symlink $file";
    #    }
    #}
}

############################################################################
# Main

# Change dir to submission directory obtained from command-line

my $wd = cwd();

GetOptions("configDir=s", "verbose=i") || usage();
$opt_verbose = 1 if (!defined $opt_verbose);
if(@ARGV != 2) {
    usage();
}

my $pipelineInstance = $ARGV[0];	# currently not used
my $submitDir = $ARGV[1];	# directory where data files are
my $configPath;
if (defined $opt_configDir) {
    if ($opt_configDir =~ /^\//) {
        $configPath = $opt_configDir;
    } else {
        $configPath = "$wd/$opt_configDir";
    }
} else {
    $configPath = "$submitDir/../config"
}


my $fields = Encode::getFields($configPath);
my $daf = Encode::getDaf($submitDir, $fields, $pipelineInstance);
my $downloadDir = Encode::downloadDir($daf);

chdir($submitDir) || die "Couldn't chdir to '$submitDir'";

my $unloadRa = 'out/unload.ra';
if(!(-e $unloadRa)) {
    HgAutomate::verbose(2, "Skipping unload b/c '$unloadRa' doesn't exist\n");
    exit(0);
}

HgAutomate::verbose(2, "Unloading project in directory $submitDir\n");

# Unload resources listed in unload.ra
my %ra = RAFile::readRaFile($unloadRa, 'tablename');
my $db;
for my $key (keys %ra) {
    my $h = $ra{$key};
    my $tablename = $h->{tablename};
    my $files = $h->{files};
    my @files = split(/\s+/, $files);

    my $str = "\nkeyword: $key\n";
    for my $field (qw(tablename type assembly files)) {
        if($h->{$field}) {
            $str .= "$field: " . $h->{$field} . "\n";
        }
    }
    $str .= "\n";
    HgAutomate::verbose(3, $str);

    my $assembly = $h->{assembly};
    if(!defined($db)) {
        $db = HgDb->new(DB => $assembly);
    }

    HgAutomate::verbose(2, "Dropping table '$tablename'\n");

    my %extendedTypes = map { $_ => 1 } @Encode::extendedTypes;
    my $type = $h->{type};
    if (exists($h->{downloadOnly}) and $h->{downloadOnly}) { # dont unload stuff which is never loaded
    } elsif($type eq "genePred" || $type =~ /^bed/ || $type eq "gtf" || $extendedTypes{$type}) {
        genericUnload($assembly, $db, $tablename);
    } elsif ($type eq "wig") {
        unloadWig($assembly, $db, $tablename);
    } elsif ($type eq "bam") {
        unloadBam($assembly, $db, $tablename);
#        unlink "$downloadDir/gbdb/$tablename.bam";
#        unlink "$downloadDir/gbdb/$tablename.bam.bai";
        unlink "$downloadDir/$tablename.bam";
        unlink "$downloadDir/$tablename.bam.bai";
    } elsif ($type eq "bigWig") {
        unloadBigWig($assembly, $db, $tablename);
#        unlink "$downloadDir/gbdb/$tablename.bw";
        unlink "$downloadDir/$tablename.bigWig";
    } elsif ($type eq "bigBed") {
        unloadBigBed($assembly, $db, $tablename);
        unlink "$downloadDir/$tablename.bigBed";
    } else {
        die "ERROR: unknown type: $h->{type} in load.ra ($PROG)\n";
    }

    # delete the download files
    my $target = $downloadDir . "/" . $h->{targetFile};

    # Just in case we are unloading an OLD unload.ra which does not name the targetFile
    if (!defined($h->{targetFile})) {
        $target = "$downloadDir/$tablename.$type.gz";

        if (@files == 1 && $files[0] =~ /^$Encode::autoCreatedPrefix/) {
            $target = "$downloadDir/raw/$tablename.$type.gz";
            if (! -d "$downloadDir/raw") {
                mkdir "$downloadDir/raw" or die "Could not create dir [$downloadDir/raw] error: [$!]\n";
                }
        }
        if ($type eq "bam") {
            $target = "$downloadDir/$tablename.$type";
        } elsif ($type eq "bigWig") {
            $target = "$downloadDir/$tablename.$type";
        } else {
            $target =~ s/ //g;  # removes space in ".bed 5.gz" for example
        }
    }

    unlink $target;
    if ($type eq "bam") {
        unlink "$target.bai";
    }
}

exit(0);
