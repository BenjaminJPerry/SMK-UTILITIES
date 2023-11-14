# 2022 Ben Perry
# MIT License
# Copyright (c) 2022 Ben Perry
# Version: 1.0
# Email: ben.perry@agresearch.co.nz

import sys
import os
import re

import pandas as pd
from Bio.Seq import translate
import argparse

# This script take output summary tables from the vASV workflow and computes
# allele frequencies, variant relative to reference, and report if these
# variants results in synonymous or non-synonymous mutations in the protein.

def readvASV(filePath):
    import pandas as pd
    vASVDF = pd.read_table(filePath, header = None,
                           names = ['sequence', 'asvCount',
                                    'sampleID', 'locusID'])
    return vASVDF

def readPlex(filePath):
    import pandas as pd
    plexDF = pd.read_csv(filePath)
    return plexDF

def mergevASVPlex(vASVDF, plexDF):
    vASVDFMerged = vASVDF.merge(plexDF, on = 'locusID', how = 'left')
    return vASVDFMerged

def locusTotals(vASVDF):
    vASVDF['libraryTotal'] = vASVDF['asvCount'].groupby(vASVDF['sampleID']).transform('sum')
    vASVDF['locusTotal'] = vASVDF['asvCount'].groupby(vASVDF['locusID']).transform('sum')
    return vASVDF

def filtervASVs(vASVDF, minCount):
    filtvASVDF = vASVDF[vASVDF['asvCount'] > minCount]
    filtvASVDF['locusFiltTotal'] = filtvASVDF['asvCount'].groupby(filtvASVDF['locusID']).transform('sum')
    return filtvASVDF

#TODO: include new functions for totalSampleReads, locusTotal, asvCount, freq

def vASVFrequencies(vASVDF):
    vASVDF['barcodeFreq'] = vASVDF['locusTotal'] / vASVDF['libraryTotal']
    vASVDF['asvFreq'] = vASVDF['asvCount'] / vASVDF['locusFiltTotal']
    return vASVDF

def correctedCDS(vASVDF):
    vASVDF['CDS'] = vASVDF.apply(
        lambda vASVDF: vASVDF['sequence'][vASVDF['ampliconRF'] - 1:] if
        vASVDF['isCDS'] else 'NaN', axis = 1)

    vASVDF['CDS'] = vASVDF.apply(
        lambda vASVDF: vASVDF['CDS'][:len(vASVDF['CDS']) - (len(vASVDF['CDS']) % 3)] if
        vASVDF['isCDS'] else 'NaN', axis = 1)
    return vASVDF

def correctedRef(vASVDF):
    vASVDF['refCDS'] = vASVDF.apply(
        lambda vASVDF:
        vASVDF['ampliconRef'][vASVDF['ampliconRF'] - 1:] if
        vASVDF['isCDS'] else 'NaN', axis = 1)
    vASVDF['refCDS'] = vASVDF.apply(
        lambda vASVDF:
        vASVDF['refCDS'][:len(vASVDF['refCDS']) - (len(vASVDF['refCDS']) % 3)] if
        vASVDF['isCDS'] else 'NaN', axis = 1)
    return vASVDF

def translateCDS(vASVDF):
    from Bio.Seq import translate
    vASVDF['AA'] = vASVDF['CDS'].apply(translate)
    vASVDF.loc[vASVDF['isCDS'] == False, 'AA'] = 'NaN'
    return vASVDF

def translateRef(vASVDF):
    from Bio.Seq import translate
    vASVDF['refAA'] = vASVDF['refCDS'].apply(translate)
    vASVDF.loc[vASVDF['isCDS'] == False, 'refAA'] = 'NaN'
    return vASVDF

def sequenceDiffs(sequence, reference, isCDS, offset, protein = False):
    diffstring = ''
    dif = '{R}-{I}-{S}'
    if protein:
        if isCDS:
            if len(sequence) != len(reference):
                return 'FRAMESHIFT'
            elif sequence == reference:
                return 'NODIFF'
            else:
                diffs = []
                for i in range(len(reference)):
                    if sequence[i] != reference[i]:
                        pos = i
                        pos = pos + (1 + offset)
                        diffs.append(
                            dif.format(R = reference[i],
                                       I = pos,
                                       S = sequence[i]))
                        diffstring = ';'.join(diffs)
                        i += 1
                    else:
                        i += 1
                return diffstring
        else:
            return 'NaN'
    else:
        if len(sequence) != len(reference):
            return 'INDEL'
        elif sequence == reference:
            return 'NODIFF'
        else:
            diffs = []
            for i in range(len(reference)):
                if sequence[i] != reference[i]:
                    pos = i
                    pos = pos + (1 + offset)
                    diffs.append(
                        dif.format(R = reference[i],
                                   I = pos,
                                   S = sequence[i]))
                    diffstring = ';'.join(diffs)
                    i += 1
                else:
                    i += 1
            return diffstring

def flagTargets(targets, AAsubs, CDS):
    import re
    if not CDS:
        return False
    elif AAsubs in ('FRAMESHIFT', 'NODIFF', 'NaN', 'INDEL'):
        return False
    elif len(AAsubs) == 0:
        return False
    else:
        targetsubs = targets.split(';')
        foundsubs = AAsubs.split(';')
        checks = []
        for target in targetsubs:
            for substitution in foundsubs:
                checks.append(bool(re.search(str(target+'-[A-Z]'), substitution)))
        return any(checks)

def get_options():
    description = """
    for each file, compute summary statistics and print out
    """
    long_description = """
    example usage:
    python summarise_vASV.py -i 6_summaries/{wildcard}.txt -m 0.1 -c 10 -o {wildcard}.summary.txt
    """
    parser = argparse.ArgumentParser(description = description,
                                     epilog = long_description,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile',
                        dest = 'infile',
                        type = str,
                        default = None,
                        required = True,
                        help = 'input summary file from vASV workflow')
    parser.add_argument('-m', '--minFreq',
                        dest = 'minFreq',
                        type = float,
                        default = 0.01,
                        help = 'minimum frequency for ASV to be reported; default: 0.01')
    parser.add_argument('-c', '--minCount',
                        dest='minCount',
                        type=int,
                        default=1,
                        help='minimum ASV count for reporting; default: 1')
    parser.add_argument('-o', '--outfile',
                        dest = 'outfile',
                        type = str,
                        required = True,
                        help = 'output csv file path')
    parser.add_argument('-p', '--plex',
                        dest = 'plexfile',
                        type = str,
                        required = True,
                        help = 'assay.csv file with 9 columns: '
                               'primerFwd,barcode,primerRev,locusID,ampliconRF,ampliconRef,isCDS,aaOffset,targetSubs')
    args = vars(parser.parse_args())
    return args

def main():
    args = get_options()
    px = readPlex(args['plexfile'])
    df = readvASV(args['infile'])
    dx = mergevASVPlex(df, px)

    locusTotals(dx)

    fdx = filtervASVs(dx, args['minCount'])

    vASVFrequencies(fdx)
    correctedCDS(fdx)
    correctedRef(fdx)
    translateCDS(fdx)
    translateRef(fdx)

    fdx['AAdiffs'] = fdx.apply(
        lambda x: sequenceDiffs(sequence = x['AA'],
                                reference = x['refAA'],
                                isCDS = x['isCDS'],
                                offset = x['aaOffset'],
                                protein = True), axis = 1)

    fdx['seqdiffs'] = fdx.apply(
        lambda x: sequenceDiffs(sequence = x['sequence'],
                                reference = x['ampliconRef'],
                                isCDS = x['isCDS'],
                                offset = 0, # For nucleotide position relative to amplicon start
                                protein = False), axis = 1)

    fdx['targetsIdentified'] = fdx.apply(
        lambda x: flagTargets(targets = x['targetSubs'],
                              AAsubs = x['AAdiffs'],
                              CDS = x['isCDS']), axis = 1)

    filtDiffs = fdx[fdx['asvFreq'] > args['minFreq']]
    filtDiffs = filtDiffs[['targetsIdentified','sampleID', 'locusID', 'seqdiffs', 'AAdiffs', 'asvFreq', 'asvCount',
                     'sequence', 'primerFwd', 'primerRev', 'barcode', 'isCDS',
                     'ampliconRF', 'ampliconRef', 'aaOffset', 'targetSubs',
                     'libraryTotal', 'locusTotal', 'locusFiltTotal', 'barcodeFreq',
                      'CDS', 'refCDS', 'AA', 'refAA']]
    filtDiffs.to_csv(args['outfile'], index = False)
    return

if __name__ == "__main__":
    main()
