#!/usr/bin/env python3
import hashlib
from os.path import basename
# Taken from http://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


md5s=[]
plinks=[ "$bed","$bim","$fam"]

for x in plinks:
    md5s.append((basename(x),md5sum(x)))




f=open("$out","w")
for (fn,md) in md5s:
    f.write("%s\\t%s\\n"%(fn,md))
f.close()

