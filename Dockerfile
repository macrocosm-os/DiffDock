# Stage 1: Build Environment Setup
FROM nvidia/cuda:11.7.1-devel-ubuntu22.04 AS builder

RUN apt-get update -y && apt-get install -y gcc wget curl git tar bzip2 unzip && rm -rf /var/lib/apt/lists/*

# Create a non-root user
ENV APPUSER="appuser"
ENV HOME=/home/$APPUSER
RUN useradd -m -u 1000 $APPUSER
USER $APPUSER
WORKDIR $HOME

ENV ENV_NAME="diffdock"
ENV DIR_NAME="DiffDock"

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
ENV PATH=$HOME/bin:$HOME/.local/bin:$PATH

# Copy and create Conda environment from environment.yml
ENV ENV_FILE_NAME=environment.yml
COPY --chown=$APPUSER:$APPUSER ./$ENV_FILE_NAME .
RUN ~/bin/micromamba env create --file $ENV_FILE_NAME && ~/bin/micromamba clean -afy --quiet

# Copy the entire repository into the container (including your API code, e.g., app.py)
COPY --chown=$APPUSER:$APPUSER . $HOME/$DIR_NAME

# Stage 2: Runtime Environment
FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

# Create non-root user and set working directory
ENV APPUSER="appuser"
ENV HOME=/home/$APPUSER
RUN useradd -m -u 1000 $APPUSER
USER $APPUSER
WORKDIR $HOME

ENV ENV_NAME="diffdock"
ENV DIR_NAME="DiffDock"

# Copy the micromamba environment and application code from the builder stage
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/micromamba $HOME/micromamba
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/bin $HOME/bin
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/$DIR_NAME $HOME/$DIR_NAME
WORKDIR $HOME/$DIR_NAME

# Set environment variables and initialize micromamba shell
ENV MAMBA_ROOT_PREFIX=$HOME/micromamba
ENV PATH=$HOME/bin:$HOME/.local/bin:$PATH
RUN micromamba shell init -s bash --root-prefix $MAMBA_ROOT_PREFIX

# Precompute series for SO(2) and SO(3) groups
RUN micromamba run -n ${ENV_NAME} python utils/precompute_series.py

# Expose original ports (if needed) and the API port
EXPOSE 7860 8501 8000

# Set the command to run the API server using Uvicorn.
# Make sure that your repository includes the `app.py` file with your FastAPI endpoints.
CMD ["sh", "-c", "micromamba run -n ${ENV_NAME} uvicorn app:app --host 0.0.0.0 --port 8000"]
