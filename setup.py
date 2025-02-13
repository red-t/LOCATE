from setuptools import setup, Extension
from Cython.Build import cythonize

def create_extension(name, sources):
    return Extension(
        name=name,
        sources=sources,
        libraries=["hts"]
    )

extensions = [
    create_extension("LOCATE.FileIO", ["LOCATE/FileIO.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"]),
    create_extension("LOCATE.Cluster", ["LOCATE/Cluster.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"]),
    create_extension("LOCATE.Assemble", ["LOCATE/Assemble.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"]),
    create_extension("LOCATE.Annotate", ["LOCATE/Annotate.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c", "LOCATE/src/anno_utils.c", "LOCATE/src/post_filter.c"]),
    create_extension("LOCATE.ParallelTaskExecutor", ["LOCATE/ParallelTaskExecutor.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c", "LOCATE/src/anno_utils.c", "LOCATE/src/post_filter.c", "LOCATE/src/ltr_utils.c"]),
]

setup(
    name="LOCATE",
    ext_modules=cythonize(extensions, language_level=3),
)