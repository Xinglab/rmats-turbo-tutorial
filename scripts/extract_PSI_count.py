"""
 # @author [Yuanyuan Wang]
 # @email [wyynju1993@gmail.com]
 # @create date 2020-05-25 04:35:20
 # @modify date 2020-05-25 04:35:20
 # @desc [description]
"""
# pylint: disable=import-error

from __future__ import print_function
import os, sys
from datetime import datetime

from class_exon import get_exon_class


def extract_PSI_COUNT(fn, od, samples):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "start: extract PSI and COUNT:", fn)
    f1 = os.path.join(od, os.path.splitext(os.path.basename(fn))[0] + "_PSI.txt")
    f2 = os.path.join(od, os.path.splitext(os.path.basename(fn))[0] + "_COUNT.txt")
    with open(fn, "r") as fin, open(f1, "w") as fo1, open(f2, "w") as fo2:

        exon, event_type = get_exon_class(fn)
        if len(samples) > 0:
            header = "\t".join(["ID"] + samples) + "\n"
            fo1.write(header)
            fo2.write(header)

        fin.readline()
        for line in fin:
            x = exon(line)
            PSI = "\t".join(str(i) for i in (x.IncLevel1 + x.IncLevel2))
            COUNT = "\t".join(
                str(i + j) for i, j in zip(x.IJC_SAMPLE_1 + x.IJC_SAMPLE_2, x.SJC_SAMPLE_1 + x.SJC_SAMPLE_2)
            )
            fo1.write(x.uniqID + "\t" + PSI + "\n")
            fo2.write(x.uniqID + "\t" + COUNT + "\n")
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "end  : extract PSI and COUNT:", fn)


def main():
    fn = sys.argv[1]
    od = sys.argv[2]
    samples = sys.argv[3:] if len(sys.argv) > 3 else []
    extract_PSI_COUNT(fn, od, samples)


if __name__ == "__main__":
    main()
