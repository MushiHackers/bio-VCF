import gzip, io, codecs, sys
from Bio.VCF import parser

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen


def thousandgenomes(**kwargs):
    '''This method enables to search through '1001 Genomes' database for VCF file corresponding to selected
    Arabidopsis Thaliana strain.

    Search is made based on provided strain name or strain Eco type / EC number.

    Method returns VCF.Reader object of selected VCF file stream.'''

    name = kwargs.get('name',None)
    ecotype = kwargs.get('ecotype',None)
    ecnumber = kwargs.get('ecnumber',None)

    if not (ecotype or name or ecnumber):
        raise Exception('You must provide strain name or strain eco type / EC number')
    if not ecotype:
        if name:
            pattern = ','+str(name)+','
            f = open("thaliana_strains.csv","r+")
            for line in f:
                if pattern in line:
                    l=line.split(',')
                    ecotype = l[0]
            if not ecotype:
                raise Exception('Provided name does not refer to any strain name in out database')
        if ecnumber:
            pattern = ','+str(ecnumber)+','
            f = open("thaliana_strains.csv","r+")
            for line in f:
                if pattern in line:
                    l = line.split(',')
                    ecotype = l[0]
            if not ecotype:
                raise Exception('Provided EC number does not refer to any strain in out database')
    if ecotype:
        tfile = "http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/"
        filename = "_snp_short_indel_with_quality_reference.vcf.gz"
        thetarfile = tfile + ecotype + filename
        page = urlopen(thetarfile)
        print (thetarfile)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            vcf = parser.Reader(contents,"r")
            return vcf
        elif sys.version < '3':
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            vcf = parser.Reader(gzip_f,"r")
            return vcf




