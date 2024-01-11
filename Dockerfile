# Dockerfile to create container with QIIME2 and q2-sidle plugin
# Push to d4straub/pipesidle:<VER>

FROM condaforge/mambaforge

LABEL authors="Daniel Straub" \
    description="Docker image containing Sidle 0.1.0-beta"

COPY environment.yml /

RUN mamba env create --file /environment.yml -p /opt/conda/envs/sidle-0.1.0-beta && \
    mamba clean --all --yes
RUN apt-get update && apt-get install -y procps

# Add conda installation dir to PATH
ENV PATH /opt/conda/envs/sidle-0.1.0-beta/bin:$PATH

# refresh qiime to make sure all is included
RUN qiime dev refresh-cache
# test commands
RUN qiime --version
RUN qiime fragment-insertion --version
RUN qiime sidle --version

# Instruct R processes to use these empty files instead of clashing with a local config
RUN touch .Rprofile
RUN touch .Renviron
