#!/usr/bin/env python
"""
Reads map file, chunk size
Returns file with chromosome chunk_start chunk_end
"""


def chunk_split(map_file, output, chunk_size, chrms='', chunk=''):
    '''
    Return: chunck files in the output folder
    '''
    data = [ [dat.split('\\t')[0], dat.split('\\t')[1]] for dat in open(map_file).readlines() ]
    datas = {}
    out = open(output, 'w')
    for dat in data:
        chrm = dat[0]
        myPos = int(dat[1])
        if chrm not in datas:
            datas[chrm] = []
        datas[chrm].append(myPos)
    data = {}
    if chrms != '':
        chrms = sorted([it for it in set(chrms.split(','))])
    else:
        chrms = sorted(datas)

    chunk_size = int(chunk_size)
    max_ = {}
    min_ = {}
    myPos = {}
    for chrm in chrms:
        max_[chrm] = max(datas[chrm])
        min_[chrm] = min(datas[chrm]) - (min(datas[chrm]) % 10) + 1
        myPos[chrm] = list(range(min_[chrm], max_[chrm], chunk_size))
    if chunk == '':
        for chrm in myPos:
            for pos in myPos[chrm]:
                start_ = pos
                end_ = start_ + chunk_size - 1
                if end_ >= max_[chrm]:
                    end_ = max_[chrm]
                out.writelines(','.join([str(chrm), str(start_), str(end_)]) + '\\n')
    else:
        chunks = chunk.split(',')
        if len(chunks) > 0:
            for chunk in chunks:
                cond = True
                chunk_ = chunk.split(':')
                chrm_ = chunk_[0]
                chunk__ = chunk_[1].split('-')
                chunk_start = int(chunk__[0])
                chunk_end = int(chunk__[1])
                for pos in list(range(chunk_start, chunk_end + chunk_size, chunk_size)):
                    if cond:
                        start_ = pos
                        end_ = start_ + chunk_size - 1
                        for chrm in myPos:
                            if str(chrm_) == str(chrm):
                                if start_ >= min_[chrm] and start_ <= max_[chrm]:
                                    if (end_ >= chunk_end) or (chunk_end - end_ + 1 <= chunk_size):
                                        out.writelines(str(chrm) + "," + str(start_) + "," + str(chunk_end) + "\\n")
                                        cond = False
                                    elif end_ <= chunk_end:
                                        out.writelines(str(chrm) + "," + str(start_) + "," + str(end_) + "\\n")
    out.close()

mapFile = "${mapFile}"
outputFile = "${chunkFile}"
chunk_size = "${chunk_size}"
chromosomes = "${chromosomes}"
chunk = "${chunk}"
chunk_split(mapFile, outputFile, chunk_size, chromosomes, chunk)
