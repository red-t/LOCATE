FROM condaforge/mambaforge:latest

WORKDIR /opt/locate

COPY environment.yml .

RUN mamba env create -f environment.yml && \
    mamba clean -afy

COPY . .

RUN mamba run -n locate pip install . --no-deps && \
    mamba run -n locate make clean

ENTRYPOINT ["mamba", "run", "-n", "locate", "locate"]
CMD ["--help"]
