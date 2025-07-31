#!/usr/bin/env python3
"""
Since chr15[10-20] is invalid VCF, this script changes it to chr15_10_20
Note: this is a quick temporary hack, the source tools need to be updated to
produce valid VCF
"""
import sys

for line in sys.stdin:
    l = line.find('[')
    if l > 0:
        r = line.find(']')
        if r > l:
            toks = line[l+1:r].split('-')
            if len(toks) == 2:
                line = f'{line[:l-1]}_{toks[0]}_{toks[1]}{line[r+1:]}'
    sys.stdout.write(line)
