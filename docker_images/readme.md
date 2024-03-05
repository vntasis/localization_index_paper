# Docker images used in the analysis

These Dockerfiles were used to construct the Docker images that include the necessary
tools to run the analysis. The Docker images can be accessed from
[DockerHub](https://hub.docker.com/r/vntasis/loc_idx_paper/tags).

- __Stan__: This is the image used with [stan-nf](https://github.com/vntasis/stan-nf) to run
  the [model](beta_inference/model.stan) that estimates beta. It includes CmdStan, and
  CmdStanR. Also, the [scripts](encode_data_analysis/convert2json.R) that create the
  input files for stan-nf used this image.
- __R__
  - bayes: This is the image with the necessary libraries to run the
    [scripts](encode_data_analysis/save_estimations.R) that read the output of stan-nf.
  - general: The rest of the scripts (the bulk of the analysis) were run using this
    image. It includes all the neccessary [R](https://www.r-project.org/) libraries
    (and more).
