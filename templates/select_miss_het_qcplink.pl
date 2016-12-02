#!/usr/bin/env perl 

 
                      

\$cut_het_high="${params.cut_het_high}";
\$cut_het_low="${params.cut_het_low}";
\$cut_miss="${params.cut_miss}";

open(MISSFILE,"$imiss");
open(HETFILE,"$het");
@all=<HETFILE>;
chomp(@all);
open(OUT,">$outfname");

\$line=0;
while(\$data = <MISSFILE>){
chomp(\$data);
  if(\$line>=1){
    chomp(\$data);
    @parts_miss=split(/\\s+/,\$data);
    \$missing=\$parts_miss[6];
    @parts_het=split(/\\s+/,\$all[\$line]);
    \$meanHet=sprintf("%.3f", (\$parts_het[5]-\$parts_het[3])/\$parts_het[5]);

    if(\$missing>\$cut_miss or \$meanHet>\$cut_het_high or \$meanHet<\$cut_het_low){
print OUT \$parts_miss[1],"\t",\$parts_miss[2],"\t",\$missing,"\t",\$meanHet,"\n";
}
}


++\$line;
}
