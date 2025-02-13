.PHONY: all build clean run

all: build

build:
	python setup.py build_ext -i

clean:
	rm -r build
	rm LOCATE/Annotate.c
	rm LOCATE/Assemble.c
	rm LOCATE/Cluster.c
	rm LOCATE/FileIO.c
	rm LOCATE/ParallelTaskExecutor.c
