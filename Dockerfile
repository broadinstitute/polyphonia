FROM quay.io/broadinstitute/viral-baseimage:0.1.20

LABEL maintainer "dpark@broadinstitute.org"

ENV \
    POLYPHONIA_PATH="/opt/polyphonia" \
    CONDA_DEFAULT_ENV="default" \
    MINICONDA_PATH="/opt/miniconda"
ENV \
    PATH="$POLYPHONIA_PATH:$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV/bin:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# initiate conda environment
RUN conda create -n $CONDA_DEFAULT_ENV python=3.7
RUN echo "source activate $CONDA_DEFAULT_ENV" >> ~/.bashrc
RUN hash -r

# install specific tools and dependencies via conda
COPY requirements-conda.txt /opt/docker
RUN /bin/bash -c "set -e; sync; conda install -y --quiet --file /opt/docker/requirements-conda.txt ; conda clean -y --all"

# install actual polyphonia scripts
COPY bin/* $POLYPHONIA_PATH/

# default bash prompt
CMD ["/bin/bash"]
