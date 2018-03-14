#!/usr/bin/perl -w

use strict;
use POSIX;
use List::Util qw[min max];

my $help=0;
my $pointsFile="-";
my $polyfile="cfht_planmonitor.txt";
my @x;
my @y;
my ($vnc, $vxc, $vyc) = (2,0,1);
my ($pnc, $pxc, $pyc) = (2,0,1);
my $nfiles=1;
my $vprefix="C";
my $vpostfix="";
my $pprefix="";
my $ppostfix="";

#Process arguments
for(my $i=0;$i<@ARGV;$i++)
{
    if($ARGV[$i]=~/^\-h$/) { $help=1; }
    elsif($ARGV[$i]=~/^-p=(.+)/) { $pointsFile=$1; }
    elsif($ARGV[$i]=~/^-v=(.+)/) { $polyfile=$1;}
    elsif($ARGV[$i]=~/^-vnc=(.+)/) { $vnc=$1; }
    elsif($ARGV[$i]=~/^-vxc=(.+)/) { $vxc=$1; }
    elsif($ARGV[$i]=~/^-vyc=(.+)/) { $vyc=$1; }
    elsif($ARGV[$i]=~/^-pnc=(.+)/) { $pnc=$1; }
    elsif($ARGV[$i]=~/^-pxc=(.+)/) { $pxc=$1; }
    elsif($ARGV[$i]=~/^-pyc=(.+)/) { $pyc=$1; }
    elsif($ARGV[$i]=~/^-vprefix=(.+)/) { $vprefix=$1; }
    elsif($ARGV[$i]=~/^-vpostfix=(.+)/) { $vpostfix=$1; }
    elsif($ARGV[$i]=~/^-pprefix=(.+)/) { $pprefix=$1; }
    elsif($ARGV[$i]=~/^-ppostfix=(.+)/) { $ppostfix=$1; }
}

if($help==1)
{
    print "Usage: ./pointInPoly.pl {options}\n";
    print "-v=<verticesFile> (cfht_planmonitor.txt default)\n";
    print "-p=<Pointsfile> (stdin default)\n";
    print "-pxc=<x column>=0\n";
    print "-pyc=<y column>=1\n";
    print "-pnc=<name column>=2\n";
    print "-vxc=<polygon x column>=0\n";
    print "-vyc=<polygon y column>=1\n";
    print "-vnc=<polygon name column>=2\n";
    print "-vprefix=<vertices name prefix>\n";
    print "-vpostfix=<vertices name postfix>\n";
    print "-pprefix=<points name prefix>\n";
    print "-ppostfix=<points name postfix>\n";
    die "";
}

my $vcolmin = max(abs($vxc),abs($vyc),abs($vnc));
my $pcolmin = max(abs($pxc),abs($pyc),abs($pnc));


#Read in the polygons
my @px;
my @py;
my @polygon;
my @vertcount=(0);


open(PF,"<$polyfile") || die "Could not open polygon file\n";
while(my $line=<PF>)
{
    chomp($line);
    my @data = split(' ',$line);

    #my $test = (@data>=$vcolmin && $line!~/^\s*#/);
    #print STDERR "$#data $test $line\n";

    if(@data>=$vcolmin && $line!~/^\s*#/)
    {
	my $x = $data[$vxc];
	my $y = $data[$vyc];
	my $n;
	if(@data==2)
	{
	    $n="";
	}
	else
	{	
	    $n = ($vnc=~/^\-{0,1}\d+$/?"$vprefix$data[$vnc]$vpostfix":"$vprefix$#vertcount$vpostfix");
	}
	push(@polygon,[$n,$x,$y]);
	$vertcount[-1]++;
    }
    elsif($vertcount[-1]>0)
    {
	push(@vertcount,0);
    }
}

warn "Vertcount = $#vertcount\n";

#Algortithm taken from http://alienryderflex.com/polygon/

#now work through the points
open(IN,"<$pointsFile") || die "Could not open points file $pointsFile\n";

my $pointcount=0;

while(my $line=<IN>)
{
    chomp($line);
    last unless $line=~/\S/;
    my @data = split(' ',$line);
    if(@data>=$pcolmin && $line!~/^\s*#/)
    {
	my $x = $data[$pxc];
	my $y = $data[$pyc];
	my $n;
	if(@data==2){
	    $n="";
	}
	else
	{	
	    $n = ($pnc=~/^\-{0,1}\d+$/?"$pprefix$data[$pnc]$ppostfix":"$pprefix$pointcount$ppostfix");
	}

	print "$n $x $y ";


	my $totvertex=0;
	for(my $poly=0;$poly<@vertcount;$poly++)
	{
	    my $oddnodes=0;
	    my $j=$vertcount[$poly]-1;
	    my $vj=$totvertex+$j;
	    #print STDERR "j vj = $j $vj\n";
	    #print STDERR "$poly $vertcount[$poly] $vj $#polygon\n";

	    for(my $i=0;$i<$vertcount[$poly]; $i++)
	    {
		my $vi=$totvertex+$i;

		#my $test1 = 0|$polygon[$vi][2]<$y;
		#my $test2 = 0|$polygon[$vj][2]>=$y;
		#my $test3 = 0|$polygon[$vj][2]<$y;
		#my $test4 = 0|$polygon[$vi][2]>=$y;
		#print STDERR "$vi $vj\n";
		#print STDERR "$n $x $y $polygon[$vi][0] $polygon[$vi][1] $polygon[$vi][2] $polygon[$vj][0] $polygon[$vj][1] $polygon[$vj][2]\n";
		#print STDERR "$test1 $test2 $test3 $test4\n";

		if(($polygon[$vi][2]<$y && $polygon[$vj][2]>=$y ||
		   $polygon[$vj][2]<$y && $polygon[$vi][2]>=$y) &&
		    ($polygon[$vi][1]<=$x || $polygon[$vj][1]<=$x))
		{

		    if($polygon[$vi][1]+($y-$polygon[$vi][2])/($polygon[$vj][2]-$polygon[$vi][2])*($polygon[$vj][1]-$polygon[$vi][1])<$x)
		    {
			$oddnodes = !$oddnodes;
		    } #end if below line
		} #end if in yrange
		$j=$i;
		$vj=$totvertex+$j;
	    } #end for i

	    if($oddnodes)
	    {
		#Point is in polygon
		my $outstr = $polygon[$totvertex][0];
		$outstr =~ s/Cccd/_ccd/g;
		    
		print "$outstr ";
	    }
	    $totvertex += $vertcount[$poly];
	} #end for polygon

	print "\n";
	$pointcount++;
    } #end if good point
} #end while read pointsfile



	
#    int   i, j=polyCorners-1 ;
#    bool  oddNodes=NO      ;
#
#    for (i=0; i<polyCorners; i++) {
#	    if (polyY[i]<y && polyY[j]>=y
#			||  polyY[j]<y && polyY[i]>=y) {
#		   if (polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x) {
#			  oddNodes=!oddNodes; }}
#	       j=i; }
#
#	  return oddNodes; }
