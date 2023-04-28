FROM mambaorg/micromamba:1.4.2-jammy
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

# RUN pip install scanpy
COPY --chown=$MAMBA_USER:$MAMBA_USER ./ /scomatic
WORKDIR /scomatic

RUN pip install -r requirements.txt

RUN Rscript r_requirements_install.R

USER root
RUN gunzip PoNs/PoN.scRNAseq.hg38.tsv.gz
RUN gunzip PoNs/PoN.scATACseq.hg38.tsv.gz 
RUN gunzip RNAediting/AllEditingSites.hg38.txt.gz

# RUN ln -sf scripts/*/* scripts/

ENV PATH="$PATH:/scomatic/scripts/"
