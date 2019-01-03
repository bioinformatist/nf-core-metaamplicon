FROM nfcore/base
LABEL authors="Yu Sun" \
      description="Docker image containing all requirements for nf-core/metaamplicon pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-metaamplicon-0.1.0/bin:$PATH
