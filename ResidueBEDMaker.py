import os, re, copy, csv
from warnings import warn
from datetime import datetime
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
# from collections import namedtuple
# For debug purposes.

dprint = lambda *args: None
#from types import SimpleNamespace

__version__ = '0.4'

class ResidueBEDMaker:
    """
    A class to create a bed file given a bunch of residue annototations. See `convert` for format of these.
    BED instance reads the 4.3 GB genome data into memory.
    One instance per BED file or `.reset` between new operations or the diagnositic will be wrong... no biggy if the data is clean.
    The `convert` goes through the genome once to match and covert the annotations give and returns a list of str (BED file entries).
    The `write_report` write the status report.
    The genomefile is a large (4GB for humans) file containing the complete chromosomes.
    For humans these are: (save the list to a file and download them with BatchEntrez)

        NC_000001.11
        NC_000002.12
        NC_000003.12
        NC_000004.12
        NC_000005.10
        NC_000006.12
        NC_000007.14
        NC_000008.11
        NC_000009.12
        NC_000010.11
        NC_000011.10
        NC_000012.12
        NC_000013.11
        NC_000014.9
        NC_000015.10
        NC_000016.10
        NC_000017.11
        NC_000018.10
        NC_000019.10
        NC_000020.11
        NC_000021.9
        NC_000022.11
        NC_000023.11
        NC_000024.10

        :param genomefile: The genome genbank.
        :type genomefile: str

    """
    _diagnostic = {'good': [], 'untranslated': [], 'missized': [], 'missing': [], 'todo': [],
                        'unverified+': [], 'unverified-': []}

    def __init__(self, genomefile: str):
        """
        :param genomefile: The genome genbank.
        :type genomefile: str
        """
        assert os.path.exists(genomefile), 'Path for genomefile is invalid.'
        assert os.path.splitext(genomefile)[1] == '.gb', 'The file genomefile should be genbank (gb)'
        self.genomefile = genomefile
        with open( self.genomefile) as fh:
            namer = lambda chromosome: re.search(" chromosome (.*?),", chromosome.description).group(1)
            self.chromosomes = {namer(chr): chr for chr in SeqIO.parse(fh, 'gb')}
        self.diagnostic = copy.deepcopy(self._diagnostic)

    def reset(self):
        """
        diagnostic is an attribute that keep track of the errors. See is_valid method.
        :return:
        """
        self.diagnostic = copy.deepcopy(self._diagnostic)
        return self

    def write_report(self, reportfile):
        with open(reportfile, 'w', newline='\n') as fh:
            w = csv.DictWriter(fh, fieldnames=self.diagnostic.keys())
            w.writeheader()
            w.writerow({k: len(self.diagnostic[k]) for k in self.diagnostic})

    def convert(self, annotations):
        """
        :param annotations: annotations of interest {gene: symbol, x: uniprot: acc, resi start, y: resi end (opt.), name: description}
        :type annotations: List[{gene: str, uniprot: int, x: int, y: int, name: str}]
        :return: bed_entries
        """
        assert len(annotations), 'No annotations to map provided'
        if any([k not in a for a in annotations for k in ('gene', 'x', 'name')]):
            raise TypeError('The annotations should be like {gene: str, x: int, y: int, name: str}.')
        genenames = {f['gene'] for f in annotations}
        dprint(genenames)
        genecheck = set()  ## canonical?
        bed_entries = []
        featbag = {g: [] for g in genenames}
        for name, chromosome in self.chromosomes.items():
            for feat in chromosome.features:
                if feat.type != 'CDS':
                    continue  # Not a gene.
                elif 'pseudo' in feat.qualifiers:
                    continue  # Not a real gene
                elif 'translation' not in feat.qualifiers:
                    continue  # Not a coding gene
                elif feat.qualifiers['gene'][0] not in genenames:
                    continue  # Not a real coding gene we care about.
                else:
                    featbag[feat.qualifiers['gene'][0]].append(feat)
        for gene in genenames:
            if len(featbag[gene]) == 0:
                self.diagnostic['missing'] = list(genenames)
            else:
                feat = sorted(featbag[gene], key=lambda f: len(f.qualifiers['translation'][0]), reverse=True)[0]
                 # we have a winner!
                for target in annotations:
                    if target['gene'] == feat.qualifiers['gene'][0]:
                        dprint('\n# '+target['gene'])
                        dprint(target)
                        #self.show_feat(feat)
                        coordinates = self.resi_to_chr(target['x'], target['y'], feat)
                        dprint(coordinates)
                        if self.is_valid(feat, target, coordinates, chromosome):
                            for i in range(0, len(coordinates), 2):
                                x = coordinates[0 + i]
                                y = coordinates[1 + i]
                                bed_entries.append(f'chr{name}\t{x}\t{y}\t{target["name"]}')
                    #it is fine.
        return bed_entries

    def is_valid(self, feat, target, coordinates, chromosome):
        """
        This method checks whether all is good.
        :param feat: the genomic gene to use for the mapping
        :type feat: BioFeature
        :param target: see annotations in convert
        :param coordinates: a list of one or more coordinates in twos
        :param chromosome: name for the bed file.
        :return: whether it is valid
        :rtype: bool
        """
        tester_v = lambda i: self._verify(coordinates[0 + i],  coordinates[1 + i], feat, chromosome.seq) !=target['x']
        tester_l = lambda i: len(self.get_feature_aa(coordinates[0 + i],  coordinates[1 + i], feat, chromosome.seq))
        if len(feat.qualifiers['translation'][0]) < target['y']:
            issue = 'missized'
        elif any([tester_v(i)  and tester_l(i) for i in range(0, len(coordinates), 2)]):
                    if feat.location.strand > 1:
                        issue = 'unverified+'
                    else:
                        issue = 'unverified-'
        else:
            self.diagnostic['good'].append(feat.qualifiers["gene"][0])
            return True
        self.diagnostic[issue].append(feat.qualifiers["gene"][0])
        warn(f'[{issue}] {target["gene"]}{target["x"]}:{target["y"]}({target["name"]}) vs. {feat.qualifiers["gene"][0]} {len(feat.qualifiers["translation"][0])}')
        return False


    def resi_to_chr(self, from_resi, to_resi, feat):
        """
        Give a residue range and the cds location, give back the postions.
        The positions are multiple of twos adn some regions may span exons.
        No testing is done. is_valid does that.
        location is a Location object from Bio.SeqRecord
        The maths makes no sense. but seems to work well with forward.
        Reverse is dodgy. So a reverse CompoundLocation has the exons in biological order
        but each exon finishes with .start (smallest index) while .end is the actual start from the translational point of view.
        the first value of a pair is the smallest coordinate. The sencond the largest. So in the case of a reverse strand,
        the from_resi is the second, not the first.
        """

        def convert(loc, offset=0):
            """
            uses from_resi and to_resi from the method.
            :param feat Biopython feat
            :param offset:
            :return:
            """
            ## get x and y, but without checking out of bounds
            if loc.strand == +1:
                x = loc.start + ((from_resi - 1) * 3 - offset) + 0
                y = loc.start + ((to_resi - 1) * 3 - offset) + 2  # aa codons go in threes, hence + 2
            elif loc.strand == -1:
                x = loc.end - ((to_resi - 1) * 3 - offset) -2  # aa codons go in threes, hence + 2
                y = loc.end - ((from_resi - 1) * 3 - offset) + 0
            else:
                raise TypeError(f'what is a {loc.strand} strand??')
            return (x, y)
        location = feat.location
        if type(location).__name__ == 'FeatureLocation':
            # override needed if missized are allowed.
            (x, y) = convert(location)
            if x < int(location.start):
                dprint('Out of bound start')
                x = location.start
            if y > int(location.end):
                dprint('Out of bound end')
                y = location.end
            return (x, y)
        elif type(location).__name__ == 'CompoundLocation':
            partials = []
            offset = 0
            parts = location.parts
            for inner_location in parts:
                (x, y) = convert(inner_location, offset)
                if inner_location.strand == 1:
                    offset += inner_location.end - inner_location.start - 1  # number of nt covered. so + 1... but changed to -1 randomly.
                else:
                    offset += inner_location.end - inner_location.start + 0
                #in reverse strand, the smallest (end) is first (labelled start) while the larger is second. The exon order is inverted.
                dprint('offest',offset)
                if x > int(inner_location.end):
                    dprint(f'{x} to the right of end: {int(inner_location.end)}')
                    #self.show_feat(feat,diff=x)
                    continue
                elif y < int(inner_location.start):
                    dprint(f'{x} to the left of start: {int(inner_location.start)}')
                    continue
                else:
                    dprint('within exon', len(inner_location) / 3)
                    if x < int(inner_location.start):
                        x = int(inner_location.start)
                        dprint('left out of bounds')
                    if y > int(inner_location.end):
                        y = int(inner_location.end)
                        dprint('right out of bounds')
                    partials.extend((x, y))
                    if inner_location.strand == 1 and y < int(inner_location.end):
                        break
                    elif inner_location.strand == -1 and x > int(inner_location.start):
                        break
            # print('here',from_resi, to_resi, partials)
            return partials
        else:
            print(location)
            raise NotImplementedError(type(location).__name__)

    @staticmethod
    def show_feat(feat, diff=0):
        """Prints the Bio feature. for debug purposes."""
        parts = [f"{int(len(p) / 3)}AA ({p.start -diff}:{p.end-diff})" for p in feat.location.parts]
        print(f'EXONS OF {feat.qualifiers["gene"][0]}: {" ".join(parts)}')

    def _verify(self, x, y, feat, seq):
        """

        :param x:
        :param y:
        :param feat:
        :param seq:
        :return:
        """
        translation = self.get_feature_aa(x, y, feat, seq)
        # the offset is one. +1 for strand +1...-1 for -1
        return feat.qualifiers['translation'][0].find(translation) + feat.location.strand

    def get_feature_aa(self, x, y, feat, seq):
        """
        Gets the translation of the genomic x and y, based on feat and genomic seq.
        :param x: smallest postion
        :param y: largest postion
        :param feat:
        :type feat: FeatureLocation
        :param seq: Bio Sequence obj.
        :return: Traslated extracted sequence. For validation.
        :rtype: str
        """
        feature_loc = FeatureLocation(x, y, feat.location.strand)
        return str(feature_loc.extract(seq).translate())

"""
Tests
=======
"""
def test_fore():
    bed = ResidueBEDMaker('chrY.gb') #human_genome.gb
    print('loaded')
    data = [{"gene": "TBL1Y", "x": 1, "uniprot": "P08048", "y": 1, "expected": 7025084, "diff": 0, "name": "sense: first residue"},
                      {"gene": "TBL1Y", "x": 20, "uniprot": "P08048", "y": 20, "expected": 7025141, "diff": 0, "name": "sense: intron spanning"},
                      {"gene": "TBL1Y", "x": 14, "uniprot": "P08048", "y": 14, "expected": 7025124, "diff": -1, "name": "sense: real 2nd codon fist exon"},
                      {"gene": "TBL1Y", "x": 40, "uniprot": "P08048", "y": 40, "expected": 7043040, "diff": -1, "name": "sense: real 2nd codon 2nd exon"},
                      {"gene": "TBL1Y", "x": 69, "uniprot": "P08048", "y": 69, "expected": 7063898, "diff": 0, "name": "sense: something.."},
                      {"gene": "TBL1Y", "x": 75, "uniprot": "P08048", "y": 75, "expected": 7063916, "diff": 0, "name": "sense: third exon."}]
    be = bed.convert(data)
    for b in be:
        print()
        xb = dict(zip(('chr', 'from', 'to', 'name'), b.split('\t')))
        matched = [d for d in data if d["name"] == xb["name"]][0]
        print(b, matched["expected"])
        if int(xb["from"]) - matched["expected"] - matched["diff"] == 0  and int(xb["to"]) - matched["expected"] - matched["diff"] == 0 + 2 :
            print('... correct.')
        elif 'intron spanning' in matched["name"]:
            print('...cannot test accuracy of exon spannng. It s right?')
        else:
            print(f'...error maybe. FROM {int(xb["from"])} - {matched["expected"]} = {int(xb["from"]) - matched["expected"]}')
            print(f'...error maybe. TO {int(xb["to"])} - {matched["expected"]} = {int(xb["to"]) - matched["expected"]}')
    print(bed.diagnostic)

def test_rev():
    bed = ResidueBEDMaker('chrY.gb') #human_genome.gb
    print('loaded')

    data = [{"gene": "UTY", "x": 1, "uniprot": "P08048", "y": 1,"expected": 13479665, "diff": 0, "name": "rev: 1st aa"},
            {"gene": "UTY", "x": 30, "uniprot": "P08048", "y": 30,"expected": 13479578, "diff": -1, "name": "Q30E 2nd pos"},
            {"gene": "UTY", "x": 90, "uniprot": "P08048", "y": 90,"expected": 13470178, "diff": 0, "name": "S90P 1nd pos"},
           {"gene": "UTY", "x": 141, "uniprot": "P08048", "y": 141, "expected": 13414747, "diff": 0, "name": "141"},
            {"gene": "UTY", "x": 903, "uniprot": "P08048", "y": 903, "expected": 13335689, "diff": -1, "name": "Leu903"},
            {"gene": "UTY", "x": 940, "uniprot": "P08048", "y": 940, "expected": 13335579, "diff": -1, "name": "Val940"},
            {"gene": "UTY", "x": 1079, "uniprot": "P08048", "y": 1079, "expected": 13323118, "diff": 0, "name": "S1079P 1nd pos"},
            {"gene": "UTY", "x": 1155, "uniprot": "P08048", "y": 1155, "expected": 13305501, "diff": 0, "name": "Arg1155"}]
    be = bed.convert(data)
    for b in be:
        print()
        xb = dict(zip(('chr', 'from', 'to', 'name'), b.split('\t')))
        matched = [d for d in data if d["name"] == xb["name"]][0]
        print(b, matched["expected"])
        # fore if int(xb["from"]) - matched["expected"] - matched["diff"] == 0  and int(xb["to"]) - matched["expected"] - matched["diff"] == 0 + 2 :
        if int(xb["from"]) - matched["expected"] - matched["diff"] - 2 == 0 and int(xb["to"]) - matched["expected"] - matched["diff"] == 0:
            print('... correct.')
        elif 'intron spanning' in matched["name"]:
            print('...cannot test accuracy of exon spannng. It s right?')
        else:
            print(f'...error maybe. FROM {int(xb["from"])} - {matched["expected"]} = {int(xb["from"]) - matched["expected"]}')
            print(f'...error maybe. TO {int(xb["to"])} - {matched["expected"]} = {int(xb["to"]) - matched["expected"]}')

    print(bed.diagnostic)

if __name__ == '__main__':
    #dprint = print
    test_fore()
    test_rev()






