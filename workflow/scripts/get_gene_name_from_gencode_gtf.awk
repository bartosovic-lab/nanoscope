BEGIN{
    FS=OFS="\t"
    } 

/^#/{
    print
    }
    
!/^#/ {
    split($9,a," "); 
    for(i=1;i<=length(a);i++) {
        if (a[i] == "gene_name")
            {print a[i+1] 
            }
        }
    }
