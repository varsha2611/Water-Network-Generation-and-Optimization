from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
	    ext_modules = [Extension("epidemic_sim", ["epidemic_sim.pyx"],
                         include_dirs=['.', 'C:\\Python27\\Lib\\site-packages\\numpy\\core\\include','/usr/local_64/epd_py25-4.0.30002-rh3-amd64/lib/python2.5/site-packages/numpy-1.1.1.0001-py2.5-linux-x86_64.egg/numpy/core/include'],
                      #debugging support
                      #extra_compile_args=["-g"],
                      #extra_link_args=["-g"],
                                )]
		)

#note: if no c libraries are used then it's possible to get better performance with
#pyximport:  http://docs.cython.org/src/userguide/tutorial.html#pyximport-cython-compilation-the-easy-way
