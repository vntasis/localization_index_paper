#!/usr/bin/awk -f

# Required input variables
# totalNM
# totalCytRatio
# Gene1NM
# Gene1CytRatio

function abs(value)
{
    return (value<0?-value:value);
}

BEGIN {
    OFS="\t"

    # Calculate total number of molecules
    totalCyt=totalNM*totalCytRatio
    totalNuc=totalNM*(1-totalCytRatio)

    # Calculate number of molecules for Gene1
    Gene1Cyt=Gene1NM*Gene1CytRatio
    Gene1Nuc=Gene1NM*(1-Gene1CytRatio)

    # Calculate expected counts for Gene1
    Counts1Who=Gene1NM/totalNM
    Counts1Cyt=Gene1Cyt/totalCyt
    Counts1Nuc=Gene1Nuc/totalNuc
    Gene1eCytRatio=Counts1Cyt/(Counts1Cyt+Counts1Nuc)

    # Calculate number of molecules for Gene2
    Gene2Cyt=totalCyt-Gene1Cyt
    Gene2Nuc=totalNuc-Gene1Nuc
    Gene2NM=Gene2Cyt+Gene2Nuc
    Gene2CytRatio=Gene2Cyt/Gene2NM

    # Calculate expected counts for Gene2
    Counts2Who=Gene2NM/totalNM
    Counts2Cyt=Gene2Cyt/totalCyt
    Counts2Nuc=Gene2Nuc/totalNuc
    Gene2eCytRatio=Counts2Cyt/(Counts2Cyt+Counts2Nuc)

    # Calculate mean absolute error
    mae=(abs(Gene1CytRatio-Gene1eCytRatio) + abs(Gene2CytRatio-Gene2eCytRatio))/2


    # Print results
    print "Gene","W","N","C","Cyt_prop"
    print "1",Gene1NM,Gene1Nuc,Gene1Cyt,Gene1CytRatio
    print "2",Gene2NM,Gene2Nuc,Gene2Cyt,Gene2CytRatio
    print "Total",totalNM,totalNuc,totalCyt,totalCytRatio
    print ""
    print "Gene","CountsW","CountsN","CountsC","Estim_Cyt_prop"
    printf "1\t%.2f\t%.2f\t%.2f\t%.2f\n",
           Counts1Who,Counts1Nuc,Counts1Cyt,Gene1eCytRatio
    printf "2\t%.2f\t%.2f\t%.2f\t%.2f\n",
           Counts2Who,Counts2Nuc,Counts2Cyt,Gene2eCytRatio
    print mae

    exit
}
