#!/usr/bin/perl -w

$USAGE="\nUsage: $0 <Infile>|<Standard input>\n

** This script splits a multi-entry GenBank file into individual files,
   each named according to the LOCUS line name.\n";


if(scalar(@ARGV)==0){
  print STDERR "<$0>\n**Reading GenBank file from STANDARD INPUT (keyboard,|,<)\n";
  print STDERR "\tIf this was unexpected, press Ctrl-d to terminate.\n";
}

if(! defined($_=<>)){ die $USAGE;}

if(! /^LOCUS/) {
  while (<>) {
    if (/^LOCUS/) {last;}
  }
}
/LOCUS\s+(\w+)/;
$gbkname = $1;
print STDERR "... now writing $gbkname\n";
open(GBKFILE,">$gbkname") || die("Cannot open $gbkname");
print GBKFILE $_;

while (<>) {
  if (/^LOCUS/) {
    close GBKFILE;
    /LOCUS\s+(\w+)/;
    $gbkname = $1;
    print STDERR "... now writing $gbkname\n";
    open(GBKFILE,">$gbkname") || die("Cannot open $gbkname");
    print GBKFILE $_;
    next;
  }
  if(/^\n$/) {next;}
  print GBKFILE $_;
}
