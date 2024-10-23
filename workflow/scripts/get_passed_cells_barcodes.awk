BEGIN {
    FS=OFS=","; 
    pass_string="\"passedMB\""
    }  
NR==1 {
    for(i=1;i<=NF;i++){
        if($i ~ pass_string){column=i}
        }
    }
$column == "TRUE" {
    print $2
    }