FROM mambaorg/micromamba:1.4.2-jammy
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libcrypt-dev \
    && rm -rf /var/lib/apt/lists/*

COPY ./ /scomatic

WORKDIR /scomatic

RUN pip install -r requirements.txt

RUN Rscript r_requirements_install.R

RUN gunzip -c PoNs/PoN.scRNAseq.hg38.tsv.gz > PoNs/PoN.scRNAseq.hg38.tsv
RUN gunzip -c PoNs/PoN.scATACseq.hg38.tsv.gz > PoNs/PoN.scATACseq.hg38.tsv
RUN gunzip -c RNAediting/AllEditingSites.hg38.txt.gz > RNAediting/AllEditingSites.hg38.txt

RUN gunzip -c PoNs/PoN.scRNAseq.hg38_nochr.tsv.gz > PoNs/PoN.scRNAseq.hg38_nochr.tsv
RUN gunzip -c PoNs/PoN.scATACseq.hg38_nochr.tsv.gz > PoNs/PoN.scATACseq.hg38_nochr.tsv
RUN gunzip -c RNAediting/AllEditingSites.hg38_nochr.txt.gz > RNAediting/AllEditingSites.hg38_nochr.txt

RUN gunzip -c example_data/chr10.fa.gz > example_data/chr10.fa
RUN samtools faidx example_data/chr10.fa

# Workaround to ensure our environment is active within singularity
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
