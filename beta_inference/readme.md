# Beta estimation

This is a [Stan](https://mc-stan.org/) model. It was used with the
[Stan-NF](https://github.com/vntasis/stan-nf) pipeline to infer _beta_, the
fraction of the total (PolyA+) RNA volume localized in the cytosol from
whole-cell, nuclear, and cytosolic RNAseq data for a sample.

As an example here is how _beta_ estimates were generated for the simulated
data:

```
nextflow -bg run vntasis/stan-nf -r 0.2 \
    --steps build-model,sample,diagnose \
    --data 'data_simulated/*.json' \
    --model model.stan \
    --outdir results_simulation_data/ \
    --chains 4 --seed 2402 \
    --numSamples 4000 --numWarmup 2000 \
    --sampleParams 'adapt delta=0.9 algorithm=hmc engine=nuts max_depth=10'
```
