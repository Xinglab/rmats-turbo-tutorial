
from __future__ import print_function
import os, sys

from class_fromGTF import get_fromGTF_class


def get_novel_ID(fn_novel):
fromGTF, event_type = get_fromGTF_class(fn_novel)
ID_novel = set()
with open(fn_novel, "r") as f:
    header = f.readline()
    for line in f:
        x = fromGTF(line)
        ID_novel.update({x.uniqID})
    return ID_novel


if __name__ == "__main__":
    # fn_novel = "./rmats/post/fromGTF.novelSpliceSite.SE.txt"
    fn_novel = sys.argv[1]

    ID_novel = get_novel_ID(fn_novel)

    for ID in ID_novel:
        print(ID)
