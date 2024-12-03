#perl! -w


open OUT, ">summary_statistics_mitodiv_simulation";
for my $i (1..10000){

    my $time = 453164; 

    print "$time\n";

    my $burnin=515000;
    my $time_w_burnin=$time+515000;
    print "$time_w_burnin\n";

    system("/home/groups/schumer/shared_bin/SLiM/slim -d TIME=$burnin -d TIME2=$time_w_burnin Simulate_mitochondrial_divergence.slim");

    my $stats=qx(perl parse_SLiM_mitochondrial_vcfs.pl population1_mito.vcf population2_mito.vcf 1 1); chomp $stats;

    print OUT "$time\t$stats\n";

    system("rm population1_mito.vcf population2_mito.vcf");

}
