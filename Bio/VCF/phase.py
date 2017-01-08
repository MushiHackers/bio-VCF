import codecs
import gzip
import sys

import re

from Bio.VCF.parser import VCFReader


class _Haplotype(object):
    """Haplotype info"""

    # TODO equals

    def __init__(self, hap):
        self.is_transmitted = None
        if len(hap.split('_')) > 1:
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
            # in the other case we do not now - stays as None

    def __str__(self):
        return "/%(name)s, transmitted: %(is_transmitted)s/" % self.__dict__


class _Sample(object):
    """Sample info"""

    # TODO equals

    def __init__(self, sample, haplotype, rsID):
        self.nucleotide = None
        self.is_unresolved = None
        self.exists = None
        self.haplotype = haplotype
        self._get_nucleotide(sample)
        self.rsID = rsID

    def _get_nucleotide(self, sample):
        if sample == '-':
            self.exists = False
            self.is_unresolved = False
        else:
            self.exists = True
            if sample.upper() not in ['A', 'T', 'G', 'C']:
                self.is_unresolved = True
            else:
                self.is_unresolved = False
            self.nucleotide = sample

    def __str__(self):
        if self.exists:
            return "%(nucleotide)s" % self.__dict__
        else:
            return "-"


class _PhasedRecord(object):
    """A set of haplotypes of particular SNP. Equivalent to a row in a Phased file.

    The standard Phased fields rsID and pos are available as properties.

    The list of genotype calls is in the ``samples`` property.

    via "ftp://ftp.ncbi.nlm.nih.gov/hapmap/phasing/2009-02_phaseIII/HapMap3_r2/readme.txt":
    - ``rsID`` contains the rsID of the SNPs were the individual genotypes have been typed.
    - ``pos`` contains the physical position of these SNPs in the particular chromosome.
    """

    # TODO equals


    def __init__(self, rsID, pos, samples=None):
        self.rsID = rsID
        self.pos = pos
        self.samples = samples or []

    def __iter__(self):
        return iter(self.samples)

    def __str__(self):
        samples_string = ''
        for sample in self.samples:
            samples_string += (str(sample) + ', ')
        samples_string = samples_string[:-2]
        return "%(rsID)s\t%(pos)s\t" % self.__dict__ + samples_string + ""


class PhasedReader(object):
    """ Reader for a phased files from HAPmap project, iterator """
    # TODO equals

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
        self._position = None  # nie wiem czy bedzie potrzebne, ale zapisze
        self.filedata = {'chrom': None, 'region': None, 'data_type': None}
        self._parse_haplotypes()
        self._parse_filedata()
        self.encoding = encoding

    def __iter__(self):
        return self

    def _parse_haplotypes(self):
        """Parse the IIDS of haplotypes stored in Phased file.

        The end users shouldn't have to use this. They can access the haplotypes
        directly with ``self.haplotypes``.
        """
        line = next(self.reader)
        if line.startswith('rsID'):
            row = self._row_pattern.split(line.rstrip())
            self._position = row[1]
            self.haplotypes = []
            for hap in row[2:]:
                self.haplotypes.append(_Haplotype(hap))

    def _parse_filedata(self):
        """Parse the information stored in the name of Phased file (if possible).

        The end users shouldn't have to use this. They can access the filename infos
        directly with ``self.filedata``.
        """
        name_string = self.filename.split('.')
        if (len(name_string) >= 6) and name_string[1:4] == ['consensus', 'qc', 'poly']:
            if name_string[4].startswith('chr') and len(name_string[4].split('_')) == 2:
                self.filedata['chrom'] = name_string[4].split('_')[0][3:]
                self.filedata['region'] = name_string[4].split('_')[1]
                if name_string[5].lower() in ['d', 'duos']:
                    self.filedata['data_type'] = 'duos'
                elif name_string[5].lower() == 'unr':
                    self.filedata['data_type'] = 'unrelated'
                elif name_string[5].lower() in ['trios', 'phased']:
                    self.filedata['data_type'] = 'trios'

    def next(self):
        """Return the next record in the file."""
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        rsID = row[0]
        pos = int(row[1])
        samples = []
        for i in range(len(row) - 2):
            samples.append(_Sample(row[2 + i], self.haplotypes[i], rsID))
        record = _PhasedRecord(rsID, pos, samples)
        return record

    __next__ = next  # Python 3.X compatibility

    def get_snp_with_specific_id(self, rsID):
        """Returns SNP with given rsID."""
        record = self.next()
        found = False
        try:
            while not found:
                if str(record.rsID) == str(rsID):
                    found = True
                    print(record)
                else:
                    record = self.next()
        except StopIteration:
            print('SNP with given rsID was not found.')
            
    def get_snp_within_range(self, pos1, pos2):
        """Returns SNP/SNPs within given range. Condition: pos1 must be smaller than pos2."""
        record = self.next()
        found = False
        try:
            while record is not None:
                if int(record.pos) >= int(pos1) and int(record.pos) < int(pos2):
                    found = True
                    print(record)
                record = self.next()
        except StopIteration:
            if not found:
                print('No SNP within given range was found.')
            else:
                pass

    def fetch(self, chrom = None, region=None, fsock= None, filename=None, compressed = None, prepend_chr=False,
                 strict_whitespace=False,  encoding = 'ascii'):
        """
        Fetches snps from VCF or from a region (positions)
        - filename is a filename of the VCF
        - region is positions in a string format 'pos1-pos2',
            ex.: '1102-49658'
        - other arguments are for a vcf reader.
        """
        if not (filename or fsock or region):
            raise Exception('You must provide at least filename or fsock or region')

        result = []

        if chrom and self.filedata['chrom'] and chrom!=self.filedata['chrom']:
            raise Exception('This file is for chrom '+str(self.filedata['chrom'])+' and you wanted to search for chrom '+chrom)

        if region:
            start,end = region.split('-')
            start = int(start)
            end = int(end)
            eof = False

            # TODO decide
            # should we treat them as unsorted as below
            # if so should we sort them first
            # or should we treat them as unsorted?

            while not eof:
                try:
                    rec = self.next()
                    if rec.pos >= start and rec.pos < end:
                        result.append(rec)
                except StopIteration:
                    eof = True
        elif filename or fsock:
            vcf = VCFReader(fsock,filename,compressed, prepend_chr,strict_whitespace,encoding)
            rec = self.next()
            vcfrec = vcf.next()
            if self.filedata['chrom']:
                pass
                # TODO tooo
                # filtrowanie chromem
                # a potem pracujemy na vcfie wyfiltrowanym chromem

            # TODO write fetch vcf

        for r in result:
            print(r)
        return result


class PhasedWriter(object):
    """Phased file writer. On Windows Python 2, open stream with 'wb'."""

    def __init__(self, stream, template):
        self.template = template
        self.stream = stream
        self._write_header()

    def _write_header(self):
        """ write header (haplotype info) """
        header = 'rsID\t' + self.template._position + '\t'
        for hap in self.template.haplotypes:
            if hap.is_transmitted:
                suffix = '_A'
            else:
                suffix = '_B'
            header += hap.name + suffix + ' '  # there is always a space at the end of the line in original files
        self.stream.write(header + '\n')

    def write_record(self, record):
        """Write a record (SNPs) to the file """
        rec = record.rsID + '\t' + str(record.pos) + '\t'
        for sample in record.samples:
            if not sample.exists:
                rec += '- '
            else:
                rec += sample.nucleotide + ' '
        self.stream.write(rec + '\n')

    def close(self):
        """Try closing the writer"""
        try:
            self.stream.close()
        except:
            pass

    def flush(self):
        """Try flushing the writer"""
        try:
            self.stream.flush()
        except:
            pass
