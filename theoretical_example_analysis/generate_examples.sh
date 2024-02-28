#!/bin/bash

set  -euo pipefail

for r in $(seq 0.2 0.05 0.8)
do
    awk -f scripts/simple_examples/simple_example_ratios.awk  \
        -v totalNM=300 -v totalCytRatio="$r" -v Gene1NM=100 -v Gene1CytRatio=0.5 | \
    tail -n1 >> mae
done

awk 'BEGIN { print "MAE" }1' mae |
    r -d -e '
        options(scipen = 999);
        pdf("mae.pdf");
        opar <- par(lwd = 3);
        mae <- X[,1];
        names(mae) <- as.character(seq(0.2,0.8,0.05));

        barplot(mae, main="", xlab="Proportion of total RNA located in Cytoplasm",
        ylab="Mean absolute error", col="white", border="black");
        dev.off();
        print(mae);'
