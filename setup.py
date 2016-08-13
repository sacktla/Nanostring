from distutils.core import setup
import py2exe
import numpy
import pandas
import sys
import glob
from bs4 import BeautifulSoup
import lxml
import lxml.etree
import bs4.builder._lxml
#import zmq.libzmq



sys.setrecursionlimit(5000)

includes = ['scipy.special._ufuncs_cxx', 'scipy.sparse.csgraph._validation',"sip", "PyQt4", "matplotlib.backends",  "matplotlib.backends.backend_qt4agg", "matplotlib.figure","pylab", "numpy","matplotlib.backends.backend_tkagg"]
excludes = ['_gtkagg', '_tkagg', '_agg2', '_cairo', '_cocoaagg', '_fltkagg', '_gtk', '_gtkcairo', ]
data_files = [(r'mpl-data', glob.glob(r'C:\Anaconda\Lib\site-packages\matplotlib\mpl-data\*.*')),
               (r'mpl-data', [r'C:\Anaconda\Lib\site-packages\matplotlib\mpl-data\matplotlibrc']),
                (r'mpl-data\images',glob.glob(r'C:\Anaconda\Lib\site-packages\matplotlib\mpl-data\images\*.*')),
                (r'mpl-data\fonts',glob.glob(r'C:\Anaconda\Lib\site-packages\matplotlib\mpl-data\fonts\*.*'))]
setup(options = {
            "py2exe":{
            "dll_excludes":["MSVCP90.dll", "HID.DLL", "w9xpopen.exe","libgdk-win32-2.0-0.dll","libgobject-2.0-0.dll"],
            "includes":includes,
            "excludes":excludes,
            "packages":["lxml"]
            }
           }, console = ['interface.py'],data_files=data_files)
           
          