FROM nfcore/base:1.9
LABEL authors="Gisela Gabernet, Simon Heumos, Alexander Peltzer" \
      description="Docker image containing all requirements for nf-core/bcellmagic pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-bcellmagic-1.3.0dev/bin:$PATH
RUN conda env export --name nf-core-bcellmagic-1.3.0dev > nf-core-bcellmagic-1.3.0dev.yml
RUN ln -s /opt/conda/envs/nf-core-bcellmagic-1.3.0dev/bin/vsearch /opt/conda/envs/nf-core-bcellmagic-1.3.0dev/bin/usearch