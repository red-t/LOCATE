.PHONY: all build clean

all: build

build:
	pip install .

clean:
	rm -r build
	rm -r LOCATE.egg-info
	rm LOCATE/Annotate.c
	rm LOCATE/Assemble.c
	rm LOCATE/Cluster.c
	rm LOCATE/FileIO.c
	rm LOCATE/ParallelTaskExecutor.c
	rm LOCATE/Main.c
