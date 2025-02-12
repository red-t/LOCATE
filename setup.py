from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext = [
    Extension(name = "LOCATE.FileIO",
              sources = ["LOCATE/FileIO.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"],
              libraries = ["hts"]),

    Extension(name = "LOCATE.Cluster",
              sources = ["LOCATE/Cluster.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"],
              libraries = ["hts"]),
    
    Extension(name = "LOCATE.Assemble",
              sources = ["LOCATE/Assemble.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c"],
              libraries = ["hts"]),
    
    Extension(name = "LOCATE.Annotate",
              sources = ["LOCATE/Annotate.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c", "LOCATE/src/anno_utils.c", "LOCATE/src/post_filter.c"],
              libraries = ["hts"]),
    
    Extension(name = "LOCATE.ParallelModule",
              sources = ["LOCATE/ParallelModule.pyx", "LOCATE/src/AIList.c", "LOCATE/src/seg_utils.c", "LOCATE/src/cluster_utils.c", "LOCATE/src/io_utils.c", "LOCATE/src/anno_utils.c", "LOCATE/src/post_filter.c", "LOCATE/src/ltr_utils.c"],
              libraries = ["hts"]),
    ]

setup(ext_modules=cythonize(ext, language_level=3))