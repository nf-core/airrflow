FROM nfcore/base:1.7
LABEL authors="Gisela Gabernet, Simon Heumos, Alexander Peltzer" \
      description="Docker image containing all requirements for nf-core/bcellmagic pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-bcellmagic-1.0.0dev/bin:$PATH
RUN ln -s /opt/conda/envs/nfcore-bcellmagic-1.0.0dev/bin/vsearch /opt/conda/envs/nfcore-bcellmagic-1.0.0dev/bin/usearch
