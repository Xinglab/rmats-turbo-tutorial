import os, sys
import numpy as np
import math

from class_exon import exon_SE, exon_RI, exon_AXSS, exon_MXE


def get_exon_class(fn):
    # *--- start get the correct exon class ----*#
    if "SE" in fn:
        exon = exon_SE
    elif "RI" in fn:
        exon = exon_RI
    elif "A3SS" in fn:
        exon = exon_AXSS
    elif "A5SS" in fn:
        exon = exon_AXSS
    elif "MXE" in fn:
        exon = exon_MXE
    # elif 'squid' in fn:
    #     exon = exon_RI_squid
    else:
        print("Wrong Type Information in Input File Name. Please Motify It.")
        sys.exit()
    # *---- end get the correct exon class ----*#
    return exon


def read_rMATS(fn, readCov, minPSI, maxPSI, sigFDR, sigDeltaPSI, bgFDR, bgWithinGroupDeltaPSI):

    exon = get_exon_class(fn)

    filtered_event_list = []
    up_event_list = []
    dn_event_list = []
    bg_event_list = []
    with open(fn, "r") as f:
        header = f.readline()
        for line in f:
            x = exon(line)

            if (
                x.averageCountSample1 >= readCov
                and x.averageCountSample2 >= readCov
                and min(x.averagePsiSample1, x.averagePsiSample2) <= maxPSI
                and max(x.averagePsiSample1, x.averagePsiSample2) >= minPSI
                and len(x.chrom) <= 5
            ):
                filtered_event_list.append(x)
                if x.FDR <= sigFDR:
                    if x.IncLevelDifference >= sigDeltaPSI:
                        dn_event_list.append(x)
                    elif x.IncLevelDifference <= -sigDeltaPSI:
                        up_event_list.append(x)
                elif (
                    x.FDR >= bgFDR
                    and max(x.IncLevel1) - min(x.IncLevel1) <= bgWithinGroupDeltaPSI
                    and max(x.IncLevel2) - min(x.IncLevel2) <= bgWithinGroupDeltaPSI
                ):
                    bg_event_list.append(x)
            # considering cases when --b2 is not provided.
            elif (
                x.averageCountSample1 >= readCov
                and math.isnan(x.averageCountSample2)
                and min(x.averagePsiSample1, x.averagePsiSample2) <= maxPSI
                and max(x.averagePsiSample1, x.averagePsiSample2) >= minPSI
                and len(x.chrom) <= 5
            ):
                filtered_event_list.append(x)

    event_dict = {
        "upregulated": up_event_list,
        "downregulated": dn_event_list,
        "background": bg_event_list,
        "filtered": filtered_event_list,
    }

    return header, event_dict


if __name__ == "__main__":
    fn_rmats = str(sys.argv[1])
    fn_up, fn_dn, fn_bg, fn_filtered = [i + os.path.basename(fn_rmats) for i in ["up_", "dn_", "bg_", "filtered_"]]

    #### FILTERING PARAMETERS ####
    # defalt filtering parameter
    readCov = 20
    minPSI = 0.05
    maxPSI = 0.95
    sigFDR = 0.01
    sigDeltaPSI = 0.05
    bgFDR = 0.5
    bgWithinGroupDeltaPSI = 0.2
    # filtering parameters from input
    if len(sys.argv) > 2:
        readCov, minPSI, maxPSI, sigFDR, sigDeltaPSI, bgFDR, bgWithinGroupDeltaPSI = [
            float(x) for x in sys.argv[2].split(",")
        ]

    #### RUN ####
    header, event_dict = read_rMATS(
        fn_rmats, readCov, minPSI, maxPSI, sigFDR, sigDeltaPSI, bgFDR, bgWithinGroupDeltaPSI
    )

    with open(fn_filtered, "w") as f:
        f.write(header)
        for event in event_dict["filtered"]:
            f.write(str(event))

    if len(event_dict["upregulated"]) > 0 or len(event_dict["downregulated"]) > 0 or len(event_dict["background"]) > 0:
        with open(fn_up, "w") as f:
            f.write(header)
            for event in event_dict["upregulated"]:
                f.write(str(event))

        with open(fn_dn, "w") as f:
            f.write(header)
            for event in event_dict["downregulated"]:
                f.write(str(event))

        with open(fn_bg, "w") as f:
            f.write(header)
            for event in event_dict["background"]:
                f.write(str(event))

