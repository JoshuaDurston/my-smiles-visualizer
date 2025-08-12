FROM continuumio/miniconda3

# Install system dependencies needed for RDKit
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Update conda and install python 3.11 and RDKit
RUN conda update -n base -c defaults conda -y && \
    conda install -c conda-forge python=3.11 rdkit=2023.3.1 -y

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

CMD ["streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]
