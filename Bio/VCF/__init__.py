#!/usr/bin/env python
"""
A VCFv4.0 and 4.1 parser for Python.

Online version of PyVCF documentation is available at http://pyvcf.rtfd.org/
"""


from VCF.parser import Reader, Writer
from VCF.parser import VCFReader, VCFWriter
from VCF.filters import Base as Filter
from VCF.parser import RESERVED_INFO, RESERVED_FORMAT
from VCF.sample_filter import SampleFilter

VERSION = '0.6.8'
