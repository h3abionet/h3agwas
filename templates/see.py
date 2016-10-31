#!/usr/bin/env python
from __future__ import print_function
import sys


if len(sys.argv) == 1:
    sys.argv=["dups.py","$inpfname","\$outfname"]

f=open("see","w")
f.write(",".join(sys.argv))

f.close()
