#!/usr/bin/env python2.7
from pandas.io import pytables
import sys

store = pytables.HDFStore(sys.argv[1])
print store
store.close()