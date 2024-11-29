BEGIN{
    FS=OFS="\t"
    } 

/^#/{
    print
    }

$3 == "gene" {
    split($9,a," "); 
    if (a[4] ~ "protein_coding" || a[4] ~ "lncRNA") { 
        if ($7 == "+") {
            $4=$4-2000
            }; 
        if ($7 =="-"){
            $5=$5+2000};
        print $0
        }
    }