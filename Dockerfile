FROM nfcore/base
LABEL authors="Gisela Gabernet, Simon Heumos, Alexander Peltzer" \
      description="Docker image containing all requirements for nf-core/bcellmagic pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get install -y zip procps ghostscript
ENV PATH /opt/conda/envs/ggabernet-bcellmagic-dev/bin:$PATH
RUN ln -s /opt/conda/envs/ggabernet-bcellmagic-dev/bin/vsearch /opt/conda/envs/ggabernet-bcellmagic-dev/bin/usearch
