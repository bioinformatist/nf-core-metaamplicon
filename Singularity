From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Yu Sun
    DESCRIPTION Singularity image containing all requirements for the nf-core/metaamplicon pipeline
    VERSION dev

%environment
    PATH=/opt/conda/envs/nf-core-metaamplicon-dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
