import codecs
import collections
import csv
import gzip
import re
import sys
import io
import codecs
import os

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

try:
    from itertools import izip, count
except ImportError:
    izip = zip
    from itertools import count

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    import pybedtools
    from pybedtools import featurefuncs as f
except ImportError:
    pybedtools = None
    f = None

try:
    import cparse
except ImportError:
    cparse = None

from Bio.VCF.model import _Call, _Record, make_calldata_tuple
from Bio.VCF.model import _Substitution, _Breakend, _SingleBreakend, _SV

# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'H3': 'Flag', 'MQ': 'Float',
    'MQ0': 'Integer', 'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag',
    'VALIDATED': 'Flag', '1000G': 'Flag',

    # Keys used for structural variants
    'IMPRECISE': 'Flag', 'NOVEL': 'Flag', 'SVTYPE': 'String',
    'SVLEN': 'Integer', 'CIPOS': 'Integer', 'CIEND': 'Integer',
    'HOMLEN': 'Integer', 'HOMSEQ': 'String', 'BKPTID': 'String',
    'MEINFO': 'String', 'METRANS': 'String', 'DGVID': 'String',
    'DBVARID': 'String', 'DBRIPID': 'String', 'MATEID': 'String',
    'PARID': 'String', 'EVENT': 'String', 'CILEN': 'Integer',
    'DPADJ': 'Integer', 'CN': 'Integer', 'CNADJ': 'Integer',
    'CICN': 'Integer', 'CICNADJ': 'Integer'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GLE': 'String', 'PL': 'Integer', 'GP': 'Float', 'GQ': 'Integer',
    'HQ': 'Integer', 'PS': 'Integer', 'PQ': 'Integer', 'EC': 'Integer',
    'MQ': 'Integer',

    # Keys used for structural variants
    'CN': 'Integer', 'CNQ': 'Float', 'CNL': 'Float', 'NQ': 'Integer',
    'HAP': 'Integer', 'AHAP': 'Integer'
}

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ['fileformat', 'fileDate', 'reference']

# Conversion between value in file and Python value
field_counts = {
    '.': None,  # Unknown number of values
    'A': -1,  # Equal to the number of alternate alleles in a given record
    'G': -2,  # Equal to the number of genotypes in a given record
    'R': -3,  # Equal to the number of alleles including reference in a given record
}

_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
_SampleInfo = collections.namedtuple('SampleInfo', ['samples', 'gt_bases', 'gt_types', 'gt_phases'])
_Contig = collections.namedtuple('Contig', ['id', 'length'])


def lf_filter(feature, location, featuretype):  # funkcja do filtrowania w fetch'u
    if ((int(feature[3]) >= location[0] and int(feature[4]) <= location[1]) and feature[2] == featuretype):
        return True
    return False


def feature_filter(feature, featuretype):
    if (feature[2] == featuretype):
        return True
    return False


def truncate_feature(feature):
    feature.chrom = "chr" + feature.chrom
    return feature

def truncate_feature2(feature):
    feature.chrom = feature.chrom[3:]
    return feature


class _vcf_metadata_parser(object):
    """Parse the metadata in the header of a VCF file."""

    def __init__(self):
        super(_vcf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AGR])?,\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.alt_pattern = re.compile(r'''\#\#ALT=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),\s*
            Number=(?P<number>-?\d+|\.|[AGR]),\s*
            Type=(?P<type>.+),\s*
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.contig_pattern = re.compile(r'''\#\#contig=<
            ID=(?P<id>[^>,]+)
            (,.*length=(?P<length>-?\d+))?
            .*
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

    def vcf_field_count(self, num_str):
        """Cast vcf header numbers to integer or None"""
        if num_str is None:
            return None
        elif num_str not in field_counts:
            # Fixed, specified number
            return int(num_str)
        else:
            return field_counts[num_str]

    def read_info(self, info_string):
        """Read a meta-information INFO line."""
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: %s" % info_string)

        num = self.vcf_field_count(match.group('number'))

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'),
                     match.group('source'), match.group('version'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        """Read a meta-information FILTER line."""
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string)

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_alt(self, alt_string):
        """Read a meta-information ALTline."""
        match = self.alt_pattern.match(alt_string)
        if not match:
            raise SyntaxError(
                "One of the ALT lines is malformed: %s" % alt_string)

        alt = _Alt(match.group('id'), match.group('desc'))

        return (match.group('id'), alt)

    def read_format(self, format_string):
        """Read a meta-information FORMAT line."""
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: %s" % format_string)

        num = self.vcf_field_count(match.group('number'))

        form = _Format(match.group('id'), num,
                       match.group('type'), match.group('desc'))

        return (match.group('id'), form)

    def read_contig(self, contig_string):
        """Read a meta-contigrmation INFO line."""
        match = self.contig_pattern.match(contig_string)
        if not match:
            raise SyntaxError(
                "One of the contig lines is malformed: %s" % contig_string)
        length = self.vcf_field_count(match.group('length'))
        contig = _Contig(match.group('id'), length)
        return (match.group('id'), contig)

    def read_meta_hash(self, meta_string):
        # assert re.match("##.+=<", meta_string)
        items = meta_string.split('=', 1)
        # Removing initial hash marks
        key = items[0].lstrip('#')
        # N.B., items can have quoted values, so cannot just split on comma
        val = OrderedDict()
        state = 0
        k = ''
        v = ''
        for c in items[1].strip('[<>]'):

            if state == 0:  # reading item key
                if c == '=':
                    state = 1  # end of key, start reading value
                else:
                    k += c  # extend key
            elif state == 1:  # reading item value
                if v == '' and c == '"':
                    v += c  # include quote mark in value
                    state = 2  # start reading quoted value
                elif c == ',':
                    val[k] = v  # store parsed item
                    state = 0  # read next key
                    k = ''
                    v = ''
                else:
                    v += c
            elif state == 2:  # reading quoted item value
                if c == '"':
                    v += c  # include quote mark in value
                    state = 1  # end quoting
                else:
                    v += c
        if k != '':
            val[k] = v
        return key, val

    def read_meta(self, meta_string):
        if re.match("##.+=<", meta_string):
            return self.read_meta_hash(meta_string)
        match = self.meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip('#'), 'none'
        return match.group('key'), match.group('val')


class Reader(object):
    """ Reader for a VCF v 4.0 file, an iterator returning ``_Record objects`` """

    def __init__(self, fsock=None, filename=None, compressed=None, prepend_chr=False,
                 strict_whitespace=False, encoding='ascii'):
        """ Create a new Reader for a VCF file.
            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``
            'prepend_chr=True' will put 'chr' before all the CHROM values, useful
            for different sources.
            'strict_whitespace=True' will split records on tabs only (as with VCF
            spec) which allows you to parse files with spaces in the sample names.
        """
        super(Reader, self).__init__()

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

        if strict_whitespace:
            self._separator = '\t'
        else:
            self._separator = '\t| +'

        self._row_pattern = re.compile(self._separator)
        self._alt_pattern = re.compile('[\[\]]')

        self.reader = (line.strip() for line in self._reader if line.strip())

        #: metadata fields from header (string or hash, depending)
        self.metadata = None
        #: INFO fields from header
        self.infos = None
        #: FILTER fields from header
        self.filters = None
        #: ALT fields from header
        self.alts = None
        #: FORMAT fields from header
        self.formats = None
        #: contig fields from header
        self.contigs = None
        self.samples = None
        self._sample_indexes = None
        self._header_lines = []
        self._column_headers = []
        self._prepend_chr = prepend_chr
        self._parse_metainfo()
        self._format_cache = {}
        self.encoding = encoding
        self._bedtool = None
        self._stream = None

    def __iter__(self):
        return self

    def _parse_metainfo(self):
        """Parse the information stored in the metainfo of the VCF.
        The end users shouldn't have to use this. They can access the metainfo
        directly with ``self.metadata``."""
        for attr in ('metadata', 'infos', 'filters', 'alts', 'contigs', 'formats'):
            setattr(self, attr, OrderedDict())

        parser = _vcf_metadata_parser()

        line = next(self.reader)
        while line.startswith('##'):
            self._header_lines.append(line)

            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self.infos[key] = val

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self.filters[key] = val

            elif line.startswith('##ALT'):
                key, val = parser.read_alt(line)
                self.alts[key] = val

            elif line.startswith('##FORMAT'):
                key, val = parser.read_format(line)
                self.formats[key] = val

            elif line.startswith('##contig'):
                key, val = parser.read_contig(line)
                self.contigs[key] = val

            else:
                key, val = parser.read_meta(line)
                if key in SINGULAR_METADATA:
                    self.metadata[key] = val
                else:
                    if key not in self.metadata:
                        self.metadata[key] = []
                    self.metadata[key].append(val)

            line = next(self.reader)

        fields = self._row_pattern.split(line[1:])
        self._column_headers = fields[:9]
        self.samples = fields[9:]
        self._sample_indexes = dict([(x, i) for (i, x) in enumerate(self.samples)])

    def _map(self, func, iterable, bad='.'):
        """``map``, but make bad values None."""
        return [func(x) if x != bad else None
                for x in iterable]

    def _parse_filter(self, filt_str):
        """Parse the FILTER field of a VCF entry into a Python list
        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        """
        if filt_str == '.':
            return None
        elif filt_str == 'PASS':
            return []
        else:
            return filt_str.split(';')

    def _parse_info(self, info_str):
        """Parse the INFO field of a VCF entry into a dictionary of Python
        types.
        """
        if info_str == '.':
            return {}

        entries = info_str.split(';')
        retdict = {}

        for entry in entries:
            entry = entry.split('=', 1)
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                try:
                    val = self._map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = self._map(float, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type in ('String', 'Character'):
                try:
                    vals = entry[1].split(',')  # commas are reserved characters indicating multiple values
                    val = self._map(str, vals)
                except IndexError:
                    entry_type = 'Flag'
                    val = True

            try:
                if self.infos[ID].num == 1 and entry_type not in ('Flag',):
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def _parse_sample_format(self, samp_fmt):
        """ Parse the format of the calls in this _Record """
        samp_fmt = make_calldata_tuple(samp_fmt.split(':'))

        for fmt in samp_fmt._fields:
            try:
                entry_type = self.formats[fmt].type
                entry_num = self.formats[fmt].num
            except KeyError:
                entry_num = None
                try:
                    entry_type = RESERVED_FORMAT[fmt]
                except KeyError:
                    entry_type = 'String'
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        return samp_fmt

    def _parse_samples(self, samples, samp_fmt, site):
        """Parse a sample entry according to the format specified in the FORMAT
        column.
        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        """

        # check whether we already know how to parse this format
        if samp_fmt not in self._format_cache:
            self._format_cache[samp_fmt] = self._parse_sample_format(samp_fmt)
        samp_fmt = self._format_cache[samp_fmt]

        if cparse:
            return cparse.parse_samples(
                self.samples, samples, samp_fmt, samp_fmt._types, samp_fmt._nums, site)

        samp_data = []
        _map = self._map

        nfields = len(samp_fmt._fields)

        for name, sample in izip(self.samples, samples):

            # parse the data for this sample
            sampdat = [None] * nfields

            for i, vals in enumerate(sample.split(':')):

                # short circuit the most common
                if samp_fmt._fields[i] == 'GT':
                    sampdat[i] = vals
                    continue
                # genotype filters are a special case
                elif samp_fmt._fields[i] == 'FT':
                    sampdat[i] = self._parse_filter(vals)
                    continue
                elif not vals or vals == ".":
                    sampdat[i] = None
                    continue

                entry_num = samp_fmt._nums[i]
                entry_type = samp_fmt._types[i]

                # we don't need to split single entries
                if entry_num == 1 or ',' not in vals:

                    if entry_type == 'Integer':
                        try:
                            sampdat[i] = int(vals)
                        except ValueError:
                            sampdat[i] = float(vals)
                    elif entry_type == 'Float':
                        sampdat[i] = float(vals)
                    else:
                        sampdat[i] = vals

                    if entry_num != 1:
                        sampdat[i] = (sampdat[i])

                    continue

                vals = vals.split(',')

                if entry_type == 'Integer':
                    try:
                        sampdat[i] = _map(int, vals)
                    except ValueError:
                        sampdat[i] = _map(float, vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdat[i] = _map(float, vals)
                else:
                    sampdat[i] = vals

            # create a call object
            call = _Call(site, name, samp_fmt(*sampdat))
            samp_data.append(call)

        return samp_data

    def _parse_alt(self, str):
        if self._alt_pattern.search(str) is not None:
            # Paired breakend
            items = self._alt_pattern.split(str)
            remoteCoords = items[1].split(':')
            chr = remoteCoords[0]
            if chr[0] == '<':
                chr = chr[1:-1]
                withinMainAssembly = False
            else:
                withinMainAssembly = True
            pos = remoteCoords[1]
            orientation = (str[0] == '[' or str[0] == ']')
            remoteOrientation = (re.search('\[', str) is not None)
            if orientation:
                connectingSequence = items[2]
            else:
                connectingSequence = items[0]
            return _Breakend(chr, pos, orientation, remoteOrientation, connectingSequence, withinMainAssembly)
        elif str[0] == '.' and len(str) > 1:
            return _SingleBreakend(True, str[1:])
        elif str[-1] == '.' and len(str) > 1:
            return _SingleBreakend(False, str[:-1])
        elif str[0] == "<" and str[-1] == ">":
            return _SV(str[1:-1])
        else:
            return _Substitution(str)

    def next(self):
        """Return the next record in the file."""
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(self._parse_alt, row[4].split(','))

        try:
            qual = int(row[5])
        except ValueError:
            try:
                qual = float(row[5])
            except ValueError:
                qual = None

        filt = self._parse_filter(row[6])
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
        else:
            if fmt == '.':
                fmt = None

        record = _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes)

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    __next__ = next  # Python 3.X compatibility


    def _vcfstream_bed(self, bed_file, verbose = True, vcf = None):
        '''Function used in fetch_bed() method in case of VCF stream.'''
        page = urlopen(self._stream)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader('utf-8')
            contents = reader(gzip_f)
            self._bedtool = pybedtools.BedTool(contents)
        else:
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            self._bedtool = pybedtools.BedTool(gzip_f)

        bed = pybedtools.BedTool(bed_file).merge()

        if bed[0].chrom[:3] != "chr" and self._bedtool[0].chrom[:3] == 'chr':
            bed = bed.each(truncate_feature)
        if bed[0].chrom[:3] == "chr" and self._bedtool[0].chrom[:3] != 'chr':
            bed = bed.each(truncate_feature2)

        features = self._bedtool.intersect(bed)
        if verbose:
            for d in features:
                print(d)

        if vcf:
            new_vcf = self.create_vcf(features, vcf)
            return new_vcf

        return features


    def fetch_bed(self, bed_file, verbose = True, vcf = None):
        """Fetches VCF file records that correspond to regions contained in a BED file.
        Fetch is based on pybedtools 'intersect' method and returns a BedTool object of selected features.
        If vcf = new_vcf_filename is provided method returns new VCF.Reader object.
        BED file must be specified and pybedtools package is required."""

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        if (self.filename == 'r' and self._stream):
            return self._vcfstream_bed(bed_file, verbose, vcf)

        if (not self.filename or (self.filename == 'r' and not self._stream)):
            raise Exception('Please provide a filename (or add stream attribute)')

        if not self._bedtool:
            self._bedtool = pybedtools.BedTool(self.filename)

        arg = False
        if self._bedtool[0].chrom[:3] != "chr":
            arg = True
        if arg:
            self._bedtool = self._bedtool.each(truncate_feature)

        bed = pybedtools.BedTool(bed_file).merge()

        arg = False
        if bed[0].chrom[:3] != "chr":
            arg = True
        if arg:
            bed = bed.each(truncate_feature)

        features = self._bedtool.intersect(bed)
        if verbose:
            for d in features:
                print(d)

        if vcf:
            new_vcf = self.create_vcf(features, vcf)
            return new_vcf

        return features

    def fetch_bed_fsock(self, stream, verbose = False, vcf=None):
        '''This fetch works exactly the same as fetch_bed(), except the BED file is not required.
        Intervals used for intersection with a VCF file are provided in stream object of chosen BED file.
        This method returns a BedTool object of selected VCF features or new VCF.Reader object
        if vcf = new_vcf_filename is specified.
        Pybedtools and stream of gzipped file are required.'''

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        if not self.filename:
            raise Exception('Please provide a filename')
        if not self._bedtool:
            self._bedtool = pybedtools.BedTool(self.filename)

        arg = False
        if self._bedtool[0].chrom[:3] != "chr":
            arg = True
        if arg:
            self._bedtool = self._bedtool.each(truncate_feature)

        thetarfile = stream
        page = urlopen(thetarfile)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            bed = pybedtools.BedTool(contents)

            arg = False
            if bed[0].chrom[:3] != "chr":
                arg = True
            if arg:
                bed = bed.each(truncate_feature)

            features = self._bedtool.intersect(bed)
            if features:
                if verbose:
                    for d in features:
                        print(d)

                if vcf:
                    new_vcf = self.create_vcf(features, vcf)
                    return new_vcf

                return features
            else:
                raise Exception("Did not find any SV in provided intervals")
        else: # sys.version < '3':
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            bed = pybedtools.BedTool(gzip_f)

            arg = False
            if bed[0].chrom[:3] != "chr":
                arg = True
            if arg:
                bed = bed.each(truncate_feature)

            features = self._bedtool.intersect(bed)
            if features:
                if verbose:
                    for d in features:
                        print(d)

                if vcf:
                    new_vcf = self.create_vcf(features, vcf)
                    return new_vcf

                return features
            else:
                raise Exception("Did not find any SV in provided intervals")

    def _vcfstream_multilocal(self, chrom, local_list, verbose = True, vcf = None):
        '''Function used in fetch_multilocal() method in case of VCF stream.'''
        page = urlopen(self._stream)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader('utf-8')
            contents = reader(gzip_f)
            self._bedtool = pybedtools.BedTool(contents)
        else:
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            self._bedtool = pybedtools.BedTool(gzip_f)

        if (self._bedtool[0].chrom[:3] != 'chr' and chrom[:3] == 'chr'):
            chrom = chrom[3:]

        if (self._bedtool[0].chrom[:3] == 'chr' and chrom[:3] != 'chr'):
            chrom = 'chr' + chrom
        local_file = 'feature_local.bed'
        g = open(local_file, "w+")
        for l in local_list:
            description = chrom + " " + str(l[0]) + " " + str(l[1])
            feature = pybedtools.BedTool(description,
                                         from_string=True)
            g.write(str(feature))
        g.seek(0)
        result = self._bedtool.intersect(local_file)
        if verbose:
            for r in result:
                print(r)

        if vcf:
            new_vcf = self.create_vcf(result, vcf)
            return new_vcf

        return result



    def fetch_multilocal(self, chrom, local_list, verbose = False, vcf = None):
        """Fetches VCF records that correspond to intervals provided in 'local_list'.
        Local_list must be a list of tuples (start, end), where start and end coordinates are in in the
        zero-based, half-open coordinate system.
        Function returns selected records as a BedTool object or new VCF.Reader object
        if vcf = vcf_new_filename is specified.
        Chrom must be specified and pybedtools package is required."""

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        if (self.filename == 'r' and self._stream):
            return self._vcfstream_multilocal(chrom, local_list, verbose, vcf)

        if (not self.filename or (self.filename == 'r' and not self._stream)):
            raise Exception('Please provide a filename (or add stream attribute)')

        if not self._bedtool:
            self._bedtool = pybedtools.BedTool(self.filename)

        if chrom[:3] != 'chr':
            chrom = 'chr' + chrom

        local_file = 'feature_local.bed'
        g = open(local_file, "w+")
        for l in local_list:
            description = chrom + " " + str(l[0]) + " " + str(l[1])
            feature = pybedtools.BedTool(description,
                                         from_string=True)
            g.write(str(feature))
        g.seek(0)
        result = self.fetch_bed(local_file)
        if verbose:
            for r in result:
                print(r)

        if vcf:
            new_vcf = self.create_vcf(result, vcf)
            return new_vcf

        return result

    def _vcfstream_fetch(self, chrom, interval = None, verbose = True, vcf = None):
        '''Function used in fetch() method in case of VCF stream.'''

        page = urlopen(self._stream)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader('utf-8')
            contents = reader(gzip_f)
            self._bedtool = pybedtools.BedTool(contents)
        else:
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            self._bedtool = pybedtools.BedTool(gzip_f)

        if (self._bedtool[0].chrom[:3] != 'chr' and chrom[:3] == 'chr'):
            chrom = chrom[3:]

        if (self._bedtool[0].chrom[:3] == 'chr' and chrom[:3] != 'chr'):
            chrom = 'chr' + chrom

        if not interval:
            end_position = 0

            for v in self._bedtool:
                if (v.chrom == chrom and int(v.end) > end_position):
                    end_position = v.end
            description = chrom + " " + str(0) + " " + str(end_position)

            feature = pybedtools.BedTool(description, from_string=True)
            page = urlopen(self._stream)
            if sys.version > '3':
                gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
                reader = codecs.getreader('utf-8')
                contents = reader(gzip_f)
                self._bedtool = pybedtools.BedTool(contents)
            else:
                gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
                self._bedtool = pybedtools.BedTool(gzip_f)

            result = self._bedtool.intersect(feature)

            if verbose:
                for r in result:
                    print(r)

            if vcf:
                new_vcf = self.create_vcf(result, vcf)
                return new_vcf

            return result

        else:
            description = chrom + " " + str(interval[0]) + " " + str(interval[1])
            feature = pybedtools.BedTool(description, from_string=True)
            result = self._bedtool.intersect(feature)
            if verbose:
                for r in result:
                    print(r)

            if vcf:
                new_vcf = self.create_vcf(result, vcf)
                return new_vcf

            return result



    def fetch(self, chrom, interval = None, verbose = True, vcf=None):
        """Fetches those records from VCF file that correspond to selected chromosome and
        fit in selected interval (if provided).
        This method creates one-line pybedtool feature based on selected chromosome (and interval = [start,stop]),
        and then uses it in pybedtools intersection method.
        Function returns BedTool object representing selected VCF records or new VCF.Reader object
        if vcf = new_vcf_filename is provided.
        Chrom must be specified and interval is optional. Pybedtools package is required."""

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        if (self.filename == 'r' and self._stream):
            return self._vcfstream_fetch(chrom, interval, verbose, vcf)

        if (not self.filename or (self.filename == 'r' and not self._stream)):
            raise Exception('Please provide a filename (or add stream attribute)')

        if not self._bedtool:
            self._bedtool = pybedtools.BedTool(self.filename)

        test = self._bedtool[0]
        if test.chrom[:3] == "chr" and chrom[:3] != 'chr':
            chrom = 'chr' + chrom

        if test.chrom[:3] != 'chr' and chrom[:3] == 'chr':
            chrom = chrom[3:]



        if not interval:
            end_position = 0

            for v in self._bedtool:
                if (v.chrom == chrom and int(v.end) > end_position):
                    end_position = v.end
            description = chrom + " " + str(0) + " " + str(end_position)
            feature = pybedtools.BedTool(description, from_string=True)


            result = self._bedtool.intersect(feature)
            if verbose:
                for r in result:
                    print(r)

            if vcf:
                new_vcf = self.create_vcf(result, vcf)
                return new_vcf

            return result

        else:
            description = chrom + " " + str(interval[0]) + " " + str(interval[1])
            feature = pybedtools.BedTool(description, from_string=True)
            result = self._bedtool.intersect(feature)
            if verbose:
                for r in result:
                    print(r)

            if vcf:
                new_vcf = self.create_vcf(result, vcf)
                return new_vcf

            return result

    def _vcfstream_gff(self, gff_file, chrom, feature_type, location = None, verbose = True, vcf = None):
        '''Function used in fetch_gff() method in case of VCF stream.'''

        gff2bedfile = 'gffbed.bed'
        page = urlopen(self._stream)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb',fileobj = page)
            reader = codecs.getreader('utf-8')
            contents = reader(gzip_f)
            self._bedtool = pybedtools.BedTool(contents)
        else:
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            self._bedtool = pybedtools.BedTool(gzip_f)

        if (self._bedtool[0].chrom[:3] != 'chr' and chrom[:3] == 'chr'):
            chrom = chrom[3:]

        if (self._bedtool[0].chrom[:3] == 'chr' and chrom[:3] != 'chr'):
            chrom = 'chr' + chrom

        is_feature = False
        if location:
            print ("Finding SV corresponding to %s and chosen position" % (feature_type))
            genes = pybedtools.BedTool(gff_file).remove_invalid().saveas()

            arg = False
            if genes[0].chrom[:3] != "chr":
                arg = True
            if arg:
                genes = genes.each(truncate_feature)

            filter_feature = genes.filter(lf_filter, location, feature_type)
            g = open(gff2bedfile, 'w+')
            for feature in filter_feature:
                feature_bed = f.gff2bed(feature)
                g.write(str(feature_bed))
                is_feature = True
            g.seek(0)
            if is_feature:
                bed = pybedtools.BedTool(gff2bedfile)
                if bed[0].chrom[:3] != "chr" and self._bedtool[0].chrom[:3] == 'chr':
                    bed = bed.each(truncate_feature)
                    features = self._bedtool.intersect(bed)
                if bed[0].chrom[:3] == 'chr' and self._bedtool[0].chrom[:3] != 'chr':
                    bed = bed.each(truncate_feature2)
                    features = self._bedtool.intersect(bed)
                else:
                    features = self._bedtool.intersect(bed)
            else:
                raise Exception('Did not find any GFF features corresponding to chosen type and/or location.')

        else:
            print ("Finding SV corresponding to %s" % (feature_type))
            genes = pybedtools.BedTool(gff_file).remove_invalid().saveas()

            arg = False
            if genes[0].chrom[:3] != "chr":
                arg = True
            if arg:
                genes = genes.each(truncate_feature)

            filter_f = genes.filter(feature_filter, feature_type)
            g = open(gff2bedfile, "w+")
            for feature in filter_f:
                feature_bed = f.gff2bed(feature)
                g.write(str(feature_bed))
                is_feature = True
            g.seek(0)
            if is_feature:
                bed = pybedtools.BedTool(gff2bedfile)
                if bed[0].chrom[:3] != "chr" and self._bedtool[0].chrom[:3] == 'chr':
                    bed = bed.each(truncate_feature)
                    features = self._bedtool.intersect(bed)
                if bed[0].chrom[:3] == 'chr' and self._bedtool[0].chrom[:3] != 'chr':
                    bed = bed.each(truncate_feature2)
                    features = self._bedtool.intersect(bed)
                else:
                    features = self._bedtool.intersect(bed)

            else:
                raise Exception('Did not find any GFF features corresponding to chosen type.')
        if verbose:
            for d in features:
                print(d)

        if vcf:
            new_vcf = self.create_vcf(features, vcf)
            return new_vcf

        return features





    def fetch_gff(self, gff_file, chrom, feature_type, location = None, verbose = False, vcf = None):
        """This method enables to select desired features from a GFF/GFF2/GFF3 file and fetch VCF records
        that correspond to position of those features. Fetch is based on pybedtools 'intersection' method and returns
        a BedTool object of chosen VCF records or new VCF.Reader object if vcf = new_vcf_filename is provided.
        Gff file, chrom and feature type are required as well as pybedtools package.
        Selection of desired features from a GFF/GFF2/GFF3 file with a specified location is possible when
        provided optional parameter 'location=[start,end]'."""

        gff2bedfile = 'gffbed.bed'

        if not pybedtools:
            raise Exception('pybedtools not available, try "pip install pybedtools"?')

        if (self.filename == 'r' and self._stream):
            return self._vcfstream_gff(gff_file, chrom, feature_type, location, verbose, vcf)

        if (not self.filename or (self.filename == 'r' and not self._stream)):
            raise Exception('Please provide a filename (or add stream attribute)')

        if not self._bedtool:
            self._bedtool = pybedtools.BedTool(self.filename)


        if chrom[:3] != 'chr':
            chrom = 'chr' + chrom

        is_feature = False
        if location:
            print ("Finding SV corresponding to %s and chosen position" % (feature_type))
            genes = pybedtools.BedTool(gff_file).remove_invalid().saveas()

            arg = False
            if genes[0].chrom[:3] != "chr":
                arg = True
            if arg:
                genes = genes.each(truncate_feature)

            filter_feature = genes.filter(lf_filter, location, feature_type)
            g = open(gff2bedfile, 'w+')
            for feature in filter_feature:
                feature_bed = f.gff2bed(feature)
                g.write(str(feature_bed))
                is_feature = True
            g.seek(0)
            if is_feature:
                features = self.fetch_bed(gff2bedfile)
            else:
                raise Exception('Did not find any GFF features corresponding to chosen type and/or location.')

        else:
            print ("Finding SV corresponding to %s" % (feature_type))
            genes = pybedtools.BedTool(gff_file).remove_invalid().saveas()

            arg = False
            if genes[0].chrom[:3] != "chr":
                arg = True
            if arg:
                genes = genes.each(truncate_feature)

            filter_f = genes.filter(feature_filter, feature_type)
            g = open(gff2bedfile, "w+")
            for feature in filter_f:
                feature_bed = f.gff2bed(feature)
                g.write(str(feature_bed))
                is_feature = True
            g.seek(0)
            if is_feature:
                features = self.fetch_bed(gff2bedfile)

            else:
                raise Exception('Did not find any GFF features corresponding to chosen type.')
        if verbose:
            for d in features:
                print(d)

        if vcf:
            new_vcf = self.create_vcf(features, vcf)
            return new_vcf

        return features

    def fetch_gff_fsock(self, stream, chrom, feature_type, location = None, verbose = False, vcf = None):
        '''This method works exactly the same as fetch_gff(), except the GFF file is not required.
        The GFF file is replaced with a stream object from chosen database.
        Method returns a BedTool object of selected VCF records or new VCF.Reader object
        if vcf = new_vcf_filename is provided.
        Pybedtools and stream of gzipped file are required.'''

        thetarfile = stream
        page = urlopen(thetarfile)
        if sys.version < '3':
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            if location:
                result = self.fetch_gff(gzip_f, chrom, feature_type, location=location)
            else:
                result = self.fetch_gff(gzip_f, chrom, feature_type)
        else: # sys.version > '3'
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            if location:
                result = self.fetch_gff(contents, chrom, feature_type, location=location)
            else:
                result = self.fetch_gff(contents, chrom, feature_type)

        if verbose:
            for r in result:
                print(r)
        if vcf:
            new_vcf = self.create_vcf(result,vcf)
            return new_vcf

        return result


    def create_vcf(self,bedtool,vcf_name):
        '''Method creates VCF.Reader object from provided pybedtools.BedTool object.
        It is based on VCF.Writer class and uses Writer.write_record() method - new VCF file is saved to
        bio-VCF/Bio/VCF directory.
        BedTool object and name for new VCF file are required.
        Method returns VCF.Reader object.'''
        vcf_reader = self
        name = vcf_name+'.vcf'
        vcf_writer = Writer(open(name,'w'),vcf_reader)
        vcf_writer.close()
        f = open(name,'a')
        for record in bedtool:
            line = ""
            for r in record.fields:
                line += str(r) + "\t"
            line = line[:-1]
            line += "\n"
            f.write(line)
        f.close()
        vcf = Reader(open(name))
        return vcf



class Writer(object):
    """VCF Writer. On Windows Python 2, open stream with 'wb'."""

    # Reverse keys and values in header field count dictionary
    counts = dict((v, k) for k, v in field_counts.items())

    def __init__(self, stream, template, lineterminator="\n"):
        self.writer = csv.writer(stream, delimiter="\t",
                                 lineterminator=lineterminator,
                                 quotechar='', quoting=csv.QUOTE_NONE)
        self.template = template
        self.stream = stream

        # Order keys for INFO fields defined in the header (undefined fields
        # get a maximum key).
        self.info_order = collections.defaultdict(
            lambda: len(template.infos),
            dict(zip(template.infos.keys(), count())))

        two = '##{key}=<ID={0},Description="{1}">\n'
        four = '##{key}=<ID={0},Number={num},Type={2},Description="{3}">\n'
        _num = self._fix_field_count
        for (key, vals) in template.metadata.items():
            if key in SINGULAR_METADATA:
                vals = [vals]
            for val in vals:
                if isinstance(val, dict):
                    values = ','.join('{0}={1}'.format(key, value)
                                      for key, value in val.items())
                    stream.write('##{0}=<{1}>\n'.format(key, values))
                else:
                    stream.write('##{0}={1}\n'.format(key, val))
        for line in template.infos.values():
            stream.write(four.format(key="INFO", *line, num=_num(line.num)))
        for line in template.formats.values():
            stream.write(four.format(key="FORMAT", *line, num=_num(line.num)))
        for line in template.filters.values():
            stream.write(two.format(key="FILTER", *line))
        for line in template.alts.values():
            stream.write(two.format(key="ALT", *line))
        for line in template.contigs.values():
            if line.length:
                stream.write('##contig=<ID={0},length={1}>\n'.format(*line))
            else:
                stream.write('##contig=<ID={0}>\n'.format(*line))

        self._write_header()

    def _write_header(self):
        # TODO: write INFO, etc
        self.stream.write('#' + '\t'.join(self.template._column_headers
                                          + self.template.samples) + '\n')

    def write_record(self, record):
        """ write a record to the file """
        ffs = self._map(str, [record.CHROM, record.POS, record.ID, record.REF]) \
              + [self._format_alt(record.ALT), record.QUAL or '.', self._format_filter(record.FILTER),
                 self._format_info(record.INFO)]
        if record.FORMAT:
            ffs.append(record.FORMAT)

        samples = [self._format_sample(record.FORMAT, sample)
                   for sample in record.samples]
        self.writer.writerow(ffs + samples)

    def flush(self):
        """Flush the writer"""
        try:
            self.stream.flush()
        except AttributeError:
            pass

    def close(self):
        """Close the writer"""
        try:
            self.stream.close()
        except AttributeError:
            pass

    def _fix_field_count(self, num_str):
        """Restore header number to original state"""
        if num_str not in self.counts:
            return num_str
        else:
            return self.counts[num_str]

    def _format_alt(self, alt):
        return ','.join(self._map(str, alt))

    def _format_filter(self, flt):
        if flt == []:
            return 'PASS'
        return self._stringify(flt, none='.', delim=';')

    def _format_info(self, info):
        if not info:
            return '.'

        def order_key(field):
            # Order by header definition first, alphabetically second.
            return self.info_order[field], field

        return ';'.join(self._stringify_pair(f, info[f]) for f in
                        sorted(info, key=order_key))

    def _format_sample(self, fmt, sample):
        if hasattr(sample.data, 'GT'):
            gt = sample.data.GT
        else:
            gt = './.' if 'GT' in fmt else ''

        result = [gt] if gt else []
        # Following the VCF spec, GT is always the first item whenever it is present.
        for field in sample.data._fields:
            value = getattr(sample.data, field)
            if field == 'GT':
                continue
            if field == 'FT':
                result.append(self._format_filter(value))
            else:
                result.append(self._stringify(value))
        return ':'.join(result)

    def _stringify(self, x, none='.', delim=','):
        if type(x) == type([]):
            return delim.join(self._map(str, x, none))
        return str(x) if x is not None else none

    def _stringify_pair(self, x, y, none='.', delim=','):
        if isinstance(y, bool):
            return str(x) if y else ""
        return "%s=%s" % (str(x), self._stringify(y, none=none, delim=delim))

    def _map(self, func, iterable, none='.'):
        """``map``, but make None values none."""
        return [func(x) if x is not None else none
                for x in iterable]


def __update_readme():
    import Bio.VCF
    file('README.rst', 'w').write(Bio.VCF.__doc__)


# backwards compatibility
VCFReader = Reader
VCFWriter = Writer
