FROM kleinstein/immcantation:2.6.0
LABEL authors="Gisela Gabernet" \
      description="Docker image containing all requirements for nf-core/bcellmagic pipeline"
RUN sudo dnf install -y bash