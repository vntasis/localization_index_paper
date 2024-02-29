FROM rocker/r-ver:4.0.4

LABEL maintainer "Vasilis F. Ntasis <vasilisntasis94@gmail.com>" \
	description "Docker image for R environment"

# Install dependencies
RUN apt-get update --fix-missing && \
	apt-get install -y \
	less \
	libz-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	libxml2-dev \
	libgsl-dev \
	cargo \
	libpng-dev \
	libjpeg-dev \
	libbz2-dev \
	liblzma-dev \
	libcairo2-dev \
	libglpk-dev \
	# X11 dependencies
	libsm-dev \
	libxt-dev \
	# Clean up
	&& apt-get autoremove \
	&& apt-get clean

# Install R packages
RUN install2.r \
	tidyverse \
	vroom \
	caret \
	mlr3 \
	gplots \
	gganimate \
	ggplotify \
	gifski \
	paletteer \
	ggsci \
	ggupset \
	ggridges \
	ggpubr \
	ggstatsplot \
	ggfortify \
	ggpointdensity \
	heatmaply \
	gt \
	devtools \
	BiocManager \
	gprofiler2 \
	feather \
	seriation \
	Rfast \
	Rtsne \
	sjstats \
	restriktor \
	idr

# Install Bioconductor packages
RUN r -e "BiocManager::install(c('edgeR','limma','DESeq2','ComplexHeatmap','clusterProfiler','Sushi','genomation','fgsea','vsn','apeglm','sva','variancePartition'))"

# Install packages from github
RUN r -e "devtools::install_github(c('stemangiola/tidyHeatmap','thomasp85/patchwork','stephens999/ashr'))"

# Add here more packages
RUN install2.r \
	tidygenomics

# Clean up downloaded packages
RUN rm -r /tmp/downloaded_packages/

# Create a docker user and set home as the working directory
RUN useradd -ms /bin/bash docker
WORKDIR /home/docker
USER docker


# Run R
CMD ["R"]
