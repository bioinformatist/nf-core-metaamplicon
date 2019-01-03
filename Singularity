From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Yu Sun
    DESCRIPTION Singularity image containing all requirements for the nf-core/metaamplicon pipeline
    VERSION 0.1.0

%environment
    PATH=/opt/conda/envs/nf-core-metaamplicon-0.1.0/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
