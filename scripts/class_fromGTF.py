
from __future__ import print_function
import math
import numpy as np
import os, sys


def get_fromGTF_class(fn):
    # *--- start get the correct exon class ----*#
    if "SE" in fn:
        fromGTF, event_type = fromGTF_SE, "SE"
    elif "RI" in fn:
        fromGTF, event_type = fromGTF_RI, "RI"
    elif "A3SS" in fn:
        fromGTF, event_type = fromGTF_AXSS, "A3SS"
    elif "A5SS" in fn:
        fromGTF, event_type = fromGTF_AXSS, "A5SS"
    elif "MXE" in fn:
        fromGTF, event_type = fromGTF_MXE, "MXE"
    # elif 'squid' in fn:
    #     exon = exon_RI_squid
    else:
        print("Wrong Type Information in Input File Name. Please Motify It.")
        sys.exit()
    # *---- end get the correct exon class ----*#
    return fromGTF, event_type


class fromGTF_SE(object):
    def __init__(self, line):
        self.line = line
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.exonStart_0base,
            self.exonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.line_list
        self.uniqID = "|".join(
            [
                self.chrom + ":" + self.exonStart_0base + "-" + self.exonEnd,
                self.strand,
                self.upstreamEE,
                self.downstreamES,
            ]
        )
        return
    
    def __str__(self):
        return str(self.line)


class fromGTF_RI(object):
    def __init__(self, line):
        self.line = line
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.riExonStart_0base,
            self.riExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.line_list
        self.uniqID = "|".join(
            [
                self.chrom + ":" + self.upstreamEE + "-" + self.downstreamES,
                self.strand,
                self.upstreamES,
                self.downstreamEE,
            ]
        )
        return
    
    def __str__(self):
        return str(self.line)


class fromGTF_AXSS(object):
    def __init__(self, line):
        self.line = line
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.longExonStart_0base,
            self.longExonEnd,
            self.shortES,
            self.shortEE,
            self.flankingES,
            self.flankingEE,
        ) = self.line_list
    
        if int(self.flankingEE) <= int(self.longExonStart_0base):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + self.flankingES + "-" + self.flankingEE,
                    self.strand,
                    self.longExonStart_0base,
                    self.shortES,
                ]
            )
        elif int(self.flankingES) >= int(self.longExonEnd):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + self.flankingES + "-" + self.flankingEE,
                    self.strand,
                    self.longExonEnd,
                    self.shortEE,
                ]
            )
        else:
            sys.exit("Error: check A5SS and A3SS file")
        return
    
    def __str__(self):
        return str(self.line)


class fromGTF_MXE(object):
    def __init__(self, line):
        self.line = line
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.stExonStart_0base,
            self.stExonEnd,
            self.ndExonStart_0base,
            self.ndExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
        ) = self.line_list

        self.uniqID = "|".join(
            [
                self.chrom
                + ":"
                + self.stExonStart_0base
                + "-"
                + self.stExonEnd
                + ":"
                + self.ndExonStart_0base
                + "-"
                + self.ndExonEnd,
                self.strand,
                self.upstreamEE,
                self.downstreamES,
            ]
        )
        return

    def __str__(self):
        return str(self.line)
