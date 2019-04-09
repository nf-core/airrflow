FROM nfcore/base
LABEL authors="simon.heumos@qbic.uni-tuebingen.de; gisela.gabernet@qbic.uni-tuebingen.de; alexander.peltzer@qbic.uni-tuebingen.de" \
      description="Docker image containing all requirements for the nfcore/bcellmagic pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-becellmagic-1.0/bin:$PATH
