from distutils.core import setup, Extension
import numpy as np

ext_modules = [ Extension('sweepMod', sources = ['sweepmodule.c']) ]

setup(
        name = 'Sweep',
        version = '0.0',
        include_dirs = [np.get_include()], #Add Include path of numpy
        ext_modules = ext_modules
      )
