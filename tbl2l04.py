#!/usr/bin/env python3
"""

Convert data table encoded in UTF-8 into data files in L04 format


"""

__author__ = "Peter Kleiweg"
__version__ = "0.1"
__date__ = "2010/03/13"

#| imports

import os
import re
import sys


def xmlescape(m):
    return '_{}_'.format(ord(m.group()))

def fname(s):
    s = re.sub(r'[^-+a-zA-Z0-9]', xmlescape, s)
    return s + '.data'

def latin1(s):
    if s:
        if s[0] == '#':
            s = '&#35;' + s[1:]
        return s.encode('iso-8859-1', 'xmlcharrefreplace')
    else:
        return b''

data = []
fp = open(sys.argv[1], 'rt', encoding='utf-8')
for line in fp:
    if not line.strip():
        continue
    if line[0] == '#':
        continue
    fnamen = line.rstrip().split('\t')[1:]
    break
labels = []
for line in fp:
    if not line.strip():
        continue
    if line[0] == '#':
        continue
    i = line.rstrip('\n').split('\t')
    labels.append(latin1(i[0]))
    data.append(i[1:])
fp.close()
os.mkdir('tmp')
os.chdir('tmp')
fp = open('LABELS.txt', 'wb')
for i in range(len(labels)):
    fp.write('{} '.format(i + 1).encode('us-ascii') + labels[i] + b'\n')
fp.close()
for i in range(len(fnamen)):
    fp = open(fname(fnamen[i]), 'wb')
    fp.write(b'%utf8\n')
    for j in range(len(labels)):
        if data[j][i]:
            fp.write(b': ' + labels[j] + b'\n- ' + data[j][i].replace(' / ', '\n- ').encode('utf-8') + b'\n')
    fp.close()
