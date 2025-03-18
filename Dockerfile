FROM rbgcsail/diffdock

# Create necessary directories if they don't exist
RUN mkdir -p /home/appuser/bin

# First, create a place to copy files temporarily
WORKDIR /tmp/diffdock_files

# Copy our files here first
COPY . /tmp/diffdock_files

# Find where the DiffDock code is located and copy our infer.py there
RUN echo '#!/bin/bash' > /tmp/find_diffdock.sh && \
    echo 'DIFFDOCK_DIR=$(find /home/appuser -name "DiffDock" -type d || find /home/appuser -name "diffdock" -type d || echo "/home/appuser/DiffDock")' >> /tmp/find_diffdock.sh && \
    echo 'echo $DIFFDOCK_DIR' >> /tmp/find_diffdock.sh && \
    chmod +x /tmp/find_diffdock.sh && \
    DIFFDOCK_DIR=$(/tmp/find_diffdock.sh) && \
    mkdir -p $DIFFDOCK_DIR && \
    cp /tmp/diffdock_files/infer.py $DIFFDOCK_DIR/ && \
    mkdir -p $DIFFDOCK_DIR/results $DIFFDOCK_DIR/test_folding && \
    echo "export DIFFDOCK_DIR=$DIFFDOCK_DIR" > /home/appuser/.diffdock_path

# Install micromamba if it's not already installed
RUN if [ ! -f /home/appuser/bin/micromamba ]; then \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /home/appuser/bin micromamba && \
    chmod +x /home/appuser/bin/micromamba; \
    fi

# Make sure environment variables are set
ENV PATH=/home/appuser/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/home/appuser/micromamba

# Initialize micromamba
RUN /home/appuser/bin/micromamba shell init -s bash --root-prefix $MAMBA_ROOT_PREFIX

# Expose ports for API and UI
EXPOSE 7860 8000 8501

# Create a simple startup script
RUN echo '#!/bin/bash' > /home/appuser/start.sh && \
    echo 'source /home/appuser/.diffdock_path' >> /home/appuser/start.sh && \
    echo 'echo ""' >> /home/appuser/start.sh && \
    echo 'echo "==============================================="' >> /home/appuser/start.sh && \
    echo 'echo "Setting up DiffDock API server..."' >> /home/appuser/start.sh && \
    echo 'echo "==============================================="' >> /home/appuser/start.sh && \
    echo 'echo ""' >> /home/appuser/start.sh && \
    echo 'cd $DIFFDOCK_DIR' >> /home/appuser/start.sh && \
    echo '/home/appuser/bin/micromamba run -n diffdock python -m uvicorn infer:app --host 0.0.0.0 --port 8000' >> /home/appuser/start.sh && \
    chmod +x /home/appuser/start.sh

# Use the startup script
CMD ["/home/appuser/start.sh"]
