#!/usr/bin/env python
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
import os
import sys

name = "ppgplot"
version = "1.1"
description = "PGPLOT interface for Python originally made by Nick Patavalis"
author = "M.A. Breddels"
author_email = "breddels@astro.rug.nl"
url = "http://www.astro.rug.nl/~breddels/"

defines = []
includes = []
extra_compile_args = []
libraries = ["cpgplot", "pgplot"]
library_dirs = []

includes.append(os.path.join(get_python_inc(plat_specific=1), "Numeric"))
if os.name == "posix":
	libraries.append("X11")
	library_dirs.append("/usr/X11R6/lib/")
	libraries.append("m")
	libraries.append("g2c")
	try:
		library_dirs.append(os.environ["PGPLOT_DIR"])
	except KeyError:
		print >>sys.stderr, "WARNING: 'PGPLOT_DIR' env variable is not defined, might not find the libcpgplot.a library"
        library_dirs.append("/usr/X11R6/lib/")

else:
	raise Exception, "os not supported"

#libraries.append(os.path.join(get_python_lib(plat_specific=1), "Numeric")+"/_numpy.so")

extension = Extension("ppgplot", ["ppgplot.c"],
	include_dirs=includes,
	libraries=libraries,
	library_dirs=library_dirs,
	define_macros=defines,
	extra_compile_args=extra_compile_args
	)



setup(name=name,
	version=version,
	description=description,
	author=author,
	author_email=author_email,
	url=url,
	ext_modules=[extension]
	)

