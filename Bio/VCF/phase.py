import codecs
import csv
import gzip
import sys

import re


class _Haplotype(object):
    """Haplotype info"""

    def __init__(self,hap):
        self.is_transmitted=None
        if len(hap.split('_'))>1:
            self.name = '_'.join(hap.split('_')[:-1])
            self._transmitted(hap.split('_')[-1])
        else:
            self.name = hap

    def _transmitted(self, suffix):
        """
        via "ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2009-02_phaseIII/HapMap3_r2/readme.txt":
        Suffixes _A and _B represent transmitted and untransmitted alleles, respectively, for a particular genotype.
        For instance, NA19009_A would be the transmitted allele of genotype NA19009 at that site.
        """
        if suffix == 'A':
            self.is_transmitted = True
        elif suffix == 'B':
            self.is_transmitted = False
            ## in the other case we do not now - stays as None

    def __str__(self):
        return "/%(name)s, transmitted: %(is_transmitted)s/" % self.__dict__

class _Sample(object):
    """Sample info"""

    def __init__(self,sample,haplotype):
        self.nucleotide = None
        self.is_unresolved = None
        self.exists = None
        self.haplotype = haplotype
        self._get_nucleotide(sample)

    def _get_nucleotide(self,sample):
        if sample == '-':
            self.exists=False
            self.is_unresolved=False
        else:
            self.exists=True
            if sample.upper() not in ['A','T','G','C']:
                self.is_unresolved=True
            else:
                self.is_unresolved=False
            self.nucleotide=sample

    def __str__(self):
        return "(%(haplotype)s : %(nucleotide)s)" % self.__dict__

class _PhasedRecord(object):
    """A set of haplotypes of particular SNP. Equivalent to a row in a Phased file.

    The standard Phased fields rsID and pos are available as properties.

    The list of genotype calls is in the ``samples`` property.

    via "ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2009-02_phaseIII/HapMap3_r2/readme.txt":
    - ``rsID`` contains the rsID of the SNPs were the individual genotypes have been typed.
    - ``pos`` contains the physical position of these SNPs in the particular chromosome.
    """
    def __init__(self, rsID, pos, samples=None):
        self.rsID = rsID
        self.pos = pos
        self.samples = samples or []

    def __iter__(self):
        return iter(self.samples)

    def __str__(self):
        samples_string = ''
        for sample in self.samples:
            samples_string += (str(sample)+', ')
        samples_string = samples_string[:-2]
        return "Record(rsID=%(rsID)s, position=%(pos)s, samples = ["+samples_string+"]" % self.__dict__



class PhasedReader(object):
    """ Reader for a phased files from HAPmap project, iterator """

    def __init__(self, filename=None, fsock=None, compressed=None, encoding='ascii'):
        """ Create a new Reader for a phased file.

            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``
        """
        super(PhasedReader, self).__init__()

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if fsock:
            self._reader = fsock
            if filename is None and hasattr(fsock, 'name'):
                filename = fsock.name
                if compressed is None:
                    compressed = filename.endswith('.gz')
        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
            self._reader = open(filename, 'rb' if compressed else 'rt')
        self.filename = filename
        if compressed:
            self._reader = gzip.GzipFile(fileobj=self._reader)
            if sys.version > '3':
                self._reader = codecs.getreader(encoding)(self._reader)

        self._row_pattern = re.compile('\t| +')

        self.reader = (line.strip() for line in self._reader if line.strip())
        self.haplotypes = None
        self._position = None #nie wiem czy bedzie potrzebne, ale zapisze
        self.filedata = {'chrom':None,'region':None,'data_type':None}
        self._parse_haplotypes()
        self._parse_filedata()
        self.encoding = encoding

    def __iter__(self):
        return self

    def _parse_haplotypes(self):
        """Parse the IIDS of haplotypes stored in Phased file.

        The end users shouldn't have to use this. They can access the haplotpes
        directly with ``self.haplotypes``.
        """
        line = next(self.reader)
        if line.startswith('rsID'):
            row = self._row_pattern.split(line.rstrip())
            self._position = row[1]
            self.haplotypes=[]
            for hap in row[2:]:
                self.haplotypes.append(_Haplotype(hap))


    def _parse_filedata(self):
        """Parse the information stored in the name of Phased file (if possible).

        The end users shouldn't have to use this. They can access the filename infos
        directly with ``self.filedata``.
        """
        name_string = self.filename.split('.')
        if (len(name_string) >= 6) and name_string[1:4]==['consensus','qc','poly']:
            if name_string[4].startswith('chr') and len(name_string[4].split('_'))==2:
                self.filedata['chrom'] = name_string[4].split('_')[0][3:]
                self.filedata['region'] = name_string[4].split('_')[1]
                if name_string[5].lower() in ['d','duos']:
                    self.filedata['data_type'] = 'duos'
                elif name_string[5].lower() == 'unr':
                    self.filedata['data_type'] = 'unrelated'
                elif name_string[5].lower() in ['trios','phased']:
                    self.filedata['data_type'] = 'trios'


    def next(self):
        """Return the next record in the file."""
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        rsID = row[0]
        pos = int(row[1])
        samples=[]
        for i in range(len(row)-2):
            samples.append(_Sample(row[2+i],self.haplotypes[i]))
        record = _PhasedRecord(rsID,pos,samples)
        return record

    __next__ = next # Python 3.X compatibility





class PhasedWriter(object):
    """Phased file writer. On Winows Python 2, open stream with 'wb'."""

    ###
    def __init__(self):
        print('not implemented yet')
    pass