#perl! -w

if(@ARGV<4){
    print "usage is: perl parse_SLiM_mitochondrial_vcfs.pl population1_mito.vcf population2_mito.vcf id_number1 id_number2\n"; exit;
}

#example:
#perl parse_SLiM_mitochondrial_vcfs.pl population1_mito.vcf population2_mito.vcf 2 4

my $vcf1=shift(@ARGV); chomp $vcf1;
open IN1, $vcf1 or die "cannot open vcf for species 1\n";

my $vcf2=shift(@ARGV); chomp $vcf2;
open IN2, $vcf2 or die "cannot open vcf for species 2\n";

my $id1=shift(@ARGV); chomp $id1;

my $id2=shift(@ARGV); chomp $id2;

my $id1_col=$id1+8;
my $id2_col=$id2+8;

#print "$id1_col\t$id2_col\n";

my $seqlength=16639;

my @pos1=(); my @data1=();
while(my $line1=<IN1>){
    chomp $line1;
    if($line1 !~ /#/){
	my @elements=split(/\t/,$line1);
	push(@pos1, $elements[1]);
	#print "$elements[1]\n";
	push(@data1, $elements[$id1_col]);
    }

}

my @pos2=(); my @data2=();
while(my $line2=<IN2>){
    chomp $line2;
    if($line2 !~ /#/){
	my @elements=split(/\t/,$line2);
	push(@pos2, $elements[1]);
	#print "$elements[1]\n";
	push(@data2, $elements[$id2_col]);
    } 
}

my $div_count=0; 
for my $i (1..$seqlength){

    my $var1=0; 
    for my $k (0..scalar(@pos1)){
	my $focal1=$pos1[$k]; chomp $focal1;
	#print "$focal1\n";
	if($i eq $focal1){
	    $var1=$data1[$k]; chomp $var1;
	}

    }

    my $var2=0;
    for my $m (0..scalar(@pos2)){
	my $focal2=$pos2[$m];

	if($i eq $focal2){
            $var2=$data2[$m]; chomp $var2;
	    #print "$i\t$var2\n";
	}

    }

    if($var1 ne $var2){
	$div_count++;
    }


}

print "$id1\t$id2\t$div_count\t$seqlength\n";
