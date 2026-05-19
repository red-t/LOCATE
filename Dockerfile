FROM condaforge/miniforge3:latest

WORKDIR /opt/locate

RUN mamba env create -n locate && \
    mamba run -n locate mamba install -c conda-forge -c bioconda -c huzr huzr::locate

ENTRYPOINT ["mamba", "run", "-n", "locate", "locate"]
CMD ["--help"]
