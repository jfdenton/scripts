#!/usr/bin/env perl

# Copyright (c) 2013, Daniel S. Standage <daniel.standage@gmail.com>
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

# This script was inspired by examples written by Eric Sayers, available from
# the Entrez Programming Utilities Help manual (see
# http://www.ncbi.nlm.nih.gov/books/NBK25498/ for examples and
# http://www.ncbi.nlm.nih.gov/books/NBK25499/ for a comprehensive reference).

use strict;
use LWP::Simple;
use Getopt::Long;

sub print_usage
{
  my $OUT = shift(@_);
  print $OUT "Usage: perl ncbi-fetch.pl [options] query
  Options:
    --batch: INT      batch size in which to download records (NBCI limits the
                      number of records that can be downloaded in a single
                      request); default is 500
    --db: STRING      database from which to fetch records; default is 'nuccore'
    -h|--help         print this help message and exit
    --out: STRING     write results to the specified file; default is terminal
                      (standard output)
    --type: STRING    the return type of the file; default is 'fasta'
    -v|--verbose      print verbose output to terminal (standard error)

  Examples:
    perl ncbi-fetch.pl --db nucleotide 'chimpanzee[orgn]+AND+biomol+mrna[prop]'
    perl ncbi-fetch.pl 'Apis mellifera[Organism] AND refseq[Filter]'
";
}

sub ncbi_search
{
  my $dbstr = shift(@_);
  my $query = shift(@_);
  
  my $url = sprintf( "%s/esearch.fcgi?db=%s&term=%s&usehistory=y",
                     "http://eutils.ncbi.nlm.nih.gov/entrez/eutils",
                     $dbstr, $query );
  my $response = get($url);
  my $webenv = $1 if($response =~ m/<WebEnv>(\S+)<\/WebEnv>/);
  my $key = $1    if($response =~ m/<QueryKey>(\d+)<\/QueryKey>/);
  my $count = $1  if($response =~ m/<Count>(\d+)<\/Count>/);
  
  unless($webenv and $key and $count)
  {
    printf(STDERR"[ncbi-fetch.pl] warning: cannot parse values from esearch\n");
    printf(STDERR"                    WebEnv: '%s'\n", $webenv);
    printf(STDERR"                  QueryKey: '%s'\n", $key);
    printf(STDERR"                     Count: '%s'\n", $count);
    die();
  }
  
  return($webenv, $key, $count);
}

sub ncbi_fetch
{
  my $dbstr  = shift(@_);
  my $webenv = shift(@_);
  my $key    = shift(@_);
  my $count  = shift(@_);
  my $batch  = shift(@_);
  my $type   = shift(@_);
  my $OUT    = shift(@_);
  my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi".
            "?db=$dbstr&WebEnv=$webenv&query_key=$key&rettype=$type"
            ."&retmode=text";
  
  for(my $i = 0; $i < $count; $i += $batch)
  {
    my $data = get($url ."&retstart=$i&retmax=$batch");
    if( $data =~ m/[Ee]rror/ or
        $data =~ m/\<ERROR\>(.+)\<\/ERROR\>/ )
    {
      printf(STDERR "[ncbi-fetch.pl] fetch error: '%s'\n", $1);
      printf(STDERR "                retrying in 10s...\n");
      sleep(10);
      $data = get($url ."&retstart=$i&retmax=$batch");
      if( $data =~ m/[Ee]rror/ or
          $data =~ m/\<ERROR\>(.+)\<\/ERROR\>/ )
      {
        printf(STDERR "[ncbi-fetch.pl] fetch error again: '%s'\n", $1);
        printf(STDERR "                giving up...\n");
        die();
      }
    }
    print $OUT $data
  }
}

# Default parameter values
my $batch   = 500;
my $dbstr   = "nuccore";
my $outfile = "";
my $OUT     = \*STDOUT;
my $type    = "fasta";
my $verbose = 0;

# Parse options, arguments
GetOptions
(
  "batch=i"   => \$batch,
  "db=s"      => \$dbstr,
  "h|help"    => sub{ print_usage(\*STDOUT); exit(0) },
  "out=s"     => \$outfile,
  "type=s"    => \$type,
  "v|verbose" => \$verbose,
);
my $query = shift(@ARGV) or do{ print_usage(\*STDERR); die() };
if($outfile ne "")
{
  open($OUT, ">", $outfile) or die("error: unable to open outfile '$outfile'");
}

# Post esearch URL
my($webenv, $key, $count) = ncbi_search($dbstr, $query);
printf(STDERR "[ncbi-fetch.pl] Found %d records matching query\n", $count);
if($verbose)
{
  printf(STDERR"[ncbi-fetch.pl] WebEnv: '%s', QueryKey: '%s'\n", $webenv, $key);
  printf(STDERR"[ncbi-fetch.pl] Query: '%s'\n", $query);
}
if($count == 0)
{
  printf(STDERR "[ncbi-fetch.pl] warning: no records found, quitting\n");
  exit(0);
}

# Download data, write to output stream
ncbi_fetch($dbstr, $webenv, $key, $count, $batch, $type, $OUT);
close($OUT);
