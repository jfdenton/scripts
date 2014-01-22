#!/usr/bin/perl

use strict;
use warnings;

my $num;
my$len;
my $min = 1000;
my $max = 0;

while (<>) {
  my @a = split /\t/;
  $num += $a[3];
  $min = $a[3] if $min > $a[3];
  $max = $a[3] if $max < $a[3];
  $len++;
}

printf "Mean coverage : %.1f\n", $num/$len;
printf "Coverage range : %d = %d\n", $min, $max;
jimmy@hahnlab-pc01:~/Desktop/biocode$ cat mean_coverage.pl 
#!/usr/bin/perl

my($num,$den)=(0,0);
while (my $cov=<STDIN>) {
    $num=$num+$cov;
    $den++;
}
$cov=$num/$den;
print "Mean Coverage = $cov\n";

