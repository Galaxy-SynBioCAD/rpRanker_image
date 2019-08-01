FROM conda/miniconda3

WORKDIR /home/

RUN apt-get --quiet update && \
    apt-get --quiet --yes dist-upgrade && \
    apt-get install --quiet --yes --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    git \
    wget \
    libxext6 \
    libxrender-dev \
    xz-utils && \
    conda update -n base -c defaults conda && \
    conda install -y -c SBMLTeam python-libsbml && \
    conda install -y -c rdkit rdkit && \
    conda install -y -c openbabel openbabel && \
    #conda update -y anaconda-navigator && \
    pip install --upgrade pip && \ 
    pip install --no-cache-dir pytest && \
    pip install --no-cache-dir cobra && \
    pip install --no-cache-dir scipy && \
    git clone --single-branch --branch master https://mdulac:towlie1988@brsforge.micalis.fr/DBT_pipeline/rpRanker.git

#Below is the code used for shared data -- remove and use rpCache if you want to use that

COPY component_contribution_data.tar.xz /home/
#TODO seperate the full input_cache with fetching the individual files from the different sources
# line MetaNetX, the git for component analysis etc...
#perhaps by creating a new file that fetches the files
#https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv
COPY input_cache.tar.xz /home/

RUN tar -xf component_contribution_data.tar.xz && \
    tar -xf input_cache.tar.xz && \
    python /home/rpRanker/rpCache.py && \
    mv cache rpRanker/ && \
    mv data rpRanker/component_contribution/ && \
    rm -rf input_cache && \
    rm input_cache.tar.xz && \
    rm component_contribution_data.tar.xz
