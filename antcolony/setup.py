from distutils.core import setup

import sys
if sys.version_info < (3,0):
    sys.exit('Sorry, Python < 3.0 is not supported')

setup(
    name        = 'antcolony',
    version     = '0.1.1', # TODO: might want to use commit ID here
    packages    = [ 'antcolony' ],
    package_dir = {
        '': '/home/zhong/ComputerScience/CS133/ant-colony'
    },
    package_data = {
        '': ['antcolony.so']
    }
)
