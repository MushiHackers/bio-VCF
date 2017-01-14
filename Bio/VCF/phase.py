import codecs
import gzip
import sys
import re
from copy import deepcopy

from Bio.VCF.parser import VCFReader

try:
    import pybedtools
except ImportError:
    pybedtools = None



class _Haplotype(object):
    """Haplotype info"""

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
            
    def __eq__(self,other):
        return self.name==self.name and self.is_transmitted == self.is_transmitted

    def __str__(self):
        return "%(name)s, transmitted: %(is_transmitted)s" % self.__dict__


class _Sample(object):
    """Sample info"""

    def __init__(self, sample, haplotype, rsID):
        self.nucleotide = None
        self.is_unresolved = None
        self.exists = None
        self.haplotype = haplotype
        self.is_not_matching_snp = None
        self._get_nucleotide(sample)
        self.rsID = rsID

    def _get_nucleotide(self, sample):
        if len(sample)>=2:
            self.is_not_matching_snp = True
            sample = sample[0]
        else:
            self.is_not_matching_snp = False

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
    def __eq__(self, other):
        return self.rsID == other.rsID and self.haplotype == other.haplotype and self.nucleotide == other.nucleotide and self.is_unresolved == other.is_unresolved and self.exists == other.exists
        # the only thing not compares it is_not_matching_snp
            
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

    def __init__(self, rsID, pos, samples=None):
        self.rsID = rsID
        self.pos = pos
        self.samples = samples or []

    def __iter__(self):
        return iter(self.samples)
    
    def __eq__(self, other):
        if len(self.samples)!= len(other.samples):
            return False
        if sorted(self.samples, key=lambda s: (s.haplotype.name, s.haplotype.is_transmitted)) != sorted(other.samples, key=lambda s: (s.haplotype.name, s.haplotype.is_transmitted)):
            return False
        return self.rsID == self.rsID and self.pos == self.pos

    def __str__(self):
        samples_string = ''
        for sample in self.samples:
            samples_string += (str(sample) + ', ')
        samples_string = samples_string[:-2]
        return "Record(%(rsID)s at %(pos)s: " % self.__dict__ + samples_string + ")"


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

        #self._row_pattern = re.compile('\t| +')

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
            row = list(filter(None, re.split("\s+",line.rstrip())))
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
        row = list(filter(None, re.split("\s+",line.rstrip())))
        rsID = row[0]
        pos = int(row[1])
        samples = []
        for i in range(len(row) - 2):
            samples.append(_Sample(row[2 + i], self.haplotypes[i], rsID))
        record = _PhasedRecord(rsID, pos, samples)
        return record

    __next__ = next  # Python 3.X compatibility

    def get_snp_with_specific_id(self, rsID):
        """Returns SNP with user-given rsID."""
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
            print('SNP with searched rsID was not found.')
            
    def get_snp_within_range(self, pos1, pos2):
        """Returns SNPs within user-given range. Condition: pos1 must be smaller than pos2."""
        record = self.next()
        found = False
        total = int()
        try:
            while record is not None:
                if int(record.pos) >= int(pos1) and int(record.pos) < int(pos2):
                    if found is False:
                        found = True
                        print('SNPs within given range:')
                    total += 1
                    print(record)
                record = self.next()
        except StopIteration:
            if not found:
                print('No SNP within given range was found.')
            else:
                print('\n' + 'In total found: ' + str(total))
                
    def get_snp_with_specific_sample(self, haplotype, nucleotide):
        """Returns SNPs containing user-given sample."""
        record = self.next()
        searched_haplotype = _Haplotype(haplotype)
        hap_found = False
        for hap in record.samples:
            if hap.haplotype.is_transmitted == searched_haplotype.is_transmitted and hap.haplotype.name == searched_haplotype.name:
                found = False
                total = int()
                hap_found = True
                hap_index = record.samples.index(hap)
                break
        if hap_found:
            print('Searched haplotype is included in the given file...')
            try:
                while record is not None:
                    if ((record.samples)[hap_index]).exists and ((record.samples)[hap_index]).is_unresolved is False:
                        if ((record.samples)[hap_index]).nucleotide == str(nucleotide):
                            if found is False:
                                found = True
                                print('SNPs with searched sample:')
                            total += 1
                            print(str(record.rsID) + '\t' + str(record.pos))
                    elif ((record.samples)[hap_index]).exists is False:
                        not_existing += 1
                    elif ((record.samples)[hap_index]).is_unresolved:
                        unresolved += 1
                    record = self.next()
            except StopIteration:
                if not found:
                    print('No SNP with searched sample was found.')
                else:
                    print('\n' + 'In total found: ' + str(total))
        else:
            print('Searched haplotype is not included in the given file.')
            
    def get_nucleotides_from_snps_in_specific_hap(self, haplotype):
        """Returns SNPs within user-given haplotype."""
        record = self.next()
        searched_haplotype = _Haplotype(haplotype)
        hap_found = False
        for hap in record.samples:
            if hap.haplotype.is_transmitted == searched_haplotype.is_transmitted and hap.haplotype.name == searched_haplotype.name:
                found = False
                total = int()
                unresolved = int()
                not_existing = int()
                hap_found = True
                hap_index = record.samples.index(hap)
                break
        if hap_found:
            print('Searched haplotype is included in the given file...')
            try:
                while record is not None:
                    if record.samples[hap_index].is_unresolved:
                        unresolved += 1
                    elif record.samples[hap_index].exists is False:
                        not_existing += 1
                    else:
                        if found is False:
                            found = True
                            print('SNPs within given haplotype:')
                        total += 1
                        print(str(record.rsID) + '\t' + str(((record.samples)[hap_index]).nucleotide))
                    record = self.next()
            except StopIteration:
                if not found:
                    print('No SNPs were found for searched haplotype.')
                    print('SNPs found unresolved: ' + str(unresolved))
                    print('SNPs found not existing: ' + str(not_existing))                       
                else:
                    print('\n' + 'In total found: ' + str(total))
                    print('\n' + 'Found unresolved: ' + str(unresolved))
                    print('Found not existing: ' + str(not_existing))
        else:
            print('Searched haplotype is not included in the given file.')
            
    def get_sample_from_specific_snp(self, rsID, haplotype):
        record = self.next()
        searched_haplotype = _Haplotype(haplotype)
        found = False
        hap_found = False
        for hap in record.samples:
            if hap.haplotype.is_transmitted == searched_haplotype.is_transmitted and hap.haplotype.name == searched_haplotype.name:
                hap_found = True
                hap_index = record.samples.index(hap)
                break
        if hap_found:
            print('Searched haplotype is included in the given file...')
            try:
                while not found:
                    if str(record.rsID) == str(rsID):
                        found = True
                        if record.samples[hap_index].is_unresolved:
                            print('The sample for searched SNP is unresolved.')
                        elif record.samples[hap_index].exists is False:
                            print('The sample for searched SNP does not exist.')
                        else:
                            print(str(((record.samples)[hap_index]).nucleotide))
                    else:
                        record = self.next()
            except StopIteration:
                print('SNP with searched rsID is not included in the given file.')
        else:
            print('Searched haplotype is not included in the given file.')

    def fetch(self, fsock=None, filename=None, compressed=None, prepend_chr=False,
              strict_whitespace=False, encoding='ascii', verbose = True, vcf = None, not_matching_snp = "."):
        """
        Fetches snps from VCF
        - filename is a filename of the VCF
        - fsock is a stream to the file
        - verbose if true prints all the fetched records
        - vcf is a filename if a  fetched working file vcf should be saved
        - not_matching_snp is a one letter string that should be put after the snp nucleotide not matching
            the alleles in vcf file.
        - other arguments are for a vcf reader.

        WE TREAT PHASED FILES AS SORTED in fetch. Keep that in mind.

        filename/fsock fetching needs pybedtools
        """
        if not (filename or fsock):
            raise Exception('You must provide at least filename or fsock')
            
        if not_matching_snp and isinstance(not_matching_snp, str) and len(not_matching_snp)>1:
            raise Exception('not_matching_snp should be string of length 1')
            
        eof=False
        result = []

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        vfile = VCFReader(fsock, filename, compressed, prepend_chr, strict_whitespace, encoding)

        if self.filedata['chrom']:
            vfile = vfile.fetch(self.filedata['chrom'], verbose = False, vcf=vcf)

        snplist=[]

        for v in vfile:
            if isinstance(vfile, pybedtools.bedtool.BedTool):
                if len(v[3]) > 1:
                    continue
                else:
                    _alt_pattern = re.compile('[\[\]]')
                    alter = v[4].split(',')
                    for alt in alter:
                        if _alt_pattern.search(alt) is not None:
                            continue
                        elif alt[0] == '.' and len(alt) > 1:
                            continue
                        elif alt[-1] == '.' and len(alt) > 1:
                            continue
                        elif alt[0] == "<" and alt[-1] == ">":
                            continue
                        elif alt not in ['A', 'C', 'G', 'T', 'N', '*']:
                            continue
                        else:
                            snplist.append(v)
            else:
                if v.is_snp:
                    snplist.append(v)

        rec = self.next()
        while not eof:
            try:
                added = False
                for v in snplist:
                    added = False
                    if isinstance(vfile, pybedtools.bedtool.BedTool):
                        start = int(v[1]) - 1
                    else:
                        start = v.start
                    if rec.pos < start or rec.pos > start+1:
                        continue
                    if rec.pos == start:
                        alleles = [v[3]] + v[4].split(',')
                        rec_string = rec.rsID + '\t' + str(rec.pos) + '\t'

                        for sample in rec.samples:
                            if not sample.exists:
                                rec_string += '- '
                            else:
                                if not_matching_snp and sample.nucleotide not in alleles:
                                    rec_string += sample.nucleotide + not_matching_snp + ' '
                                else:
                                    rec_string += sample.nucleotide + ' '
                        result.append(rec_string)
                        rec = self.next()
                        added = True
                if not added:
                    rec = self.next()
            except StopIteration:
                eof = True
        
        if verbose:
            for r in result:
                print(r)
        origin = self._reader
        self.reader, self._reader = None,None
        resultreader = deepcopy(self)
        origin.seek(0)
        self._reader = origin
        self.reader = (line.strip() for line in self._reader if line.strip())
        next(self.reader)
        resultreader.reader = (line for line in result)
        return resultreader


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

    def write_record(self, record, not_matching_snp='.'):
        """Write a record (SNPs) to the file """
        rec = record.rsID + '\t' + str(record.pos) + '\t'
        for sample in record.samples:
            if not sample.exists:
                rec += '- '
            else:
                if not_matching_snp and sample.is_not_matching_snp:
                    rec += sample.nucleotide + not_matching_snp + ' '
                else:
                    rec += sample.nucleotide + ' '
        self.stream.write(rec + '\n')

    def close(self):
        """Try closing the writer"""
        try:
            self.stream.close()
        except AttributeError:
            pass

    def flush(self):
        """Try flushing the writer"""
        try:
            self.stream.flush()
        except AttributeError:
            pass
