#!/usr/bin/bash

flux --threads 2 -x -l -s -p whole_cell_flux.par &> flux.log
