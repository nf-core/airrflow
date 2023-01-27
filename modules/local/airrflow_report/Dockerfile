FROM condaforge/mambaforge:22.9.0-3
LABEL authors="gisela.gabernet@gmail.com" \
    description="Docker image containing base requirements for the nfcore/airrflowreport"

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get install -y procps

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-airrflow-report-3.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name nf-core-airrflow-report-3.0dev > nf-core-airrflow-report-3.0dev.yml
