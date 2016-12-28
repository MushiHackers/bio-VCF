import gzip, io, codecs, sys
from Bio.VCF import parser

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
    from urllib.request import urlretrieve
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
    from urllib import urlretrieve



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
        #print (thetarfile)
        if sys.version > '3':
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            vcf = parser.Reader(contents,"r")
            vcf._stream = thetarfile
            return vcf
        elif sys.version < '3':
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            vcf = parser.Reader(gzip_f,"r")
            vcf._stream = thetarfile
            return vcf


def thousandgenomes_country(name):
    '''This method enables to search through '1001 Genomes' database for VCF files corresponding to selected origin country of
    Arabidopsis Thaliana data.

    Country name is required and must be chosen from: "USA", "FRA", "CZE", "AUT", "KGZ", "TJK", "SWE", "UK", "GER", "KAZ",
    "BEL", "CPV", "ESP", "RUS", "NED", "FIN", "SUI", "ITA", "IRL", "POR", "EST", "DEN", "IND", "LTU", "JPN", "POL", "NOR",
    "CAN", "UKR", "AZE", "GEO", "ARM", "MAR", "CRO", "BUL", "GRC", "SVK", "ROU", "UZB", "SRB", "CHN", "IRN", "LBN", "MAR",
    "AFG".

     Method returns list of VCF.Reader objectc of selected VCF file streams.'''

    f = open("thaliana_strains.csv","r+")
    f.readline()
    pattern = ','+str(name) + ','
    e_list=[]
    for line in f:
        if pattern in line:
            l = line.split(',')
            e_list.append(l[0])
    tfile = "http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/"
    filename = "_snp_short_indel_with_quality_reference.vcf.gz"
    if sys.version > '3':
        resulting_vcf = []
        for e in e_list:
            thetarfile = tfile + e + filename
            page = urlopen(thetarfile)
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            vcf = parser.Reader(contents,'r')
            vcf._stream = thetarfile
            resulting_vcf.append(vcf)
            print (resulting_vcf)
        return resulting_vcf
    elif sys.version < '3':
        resulting_vcf=[]
        for e in e_list:
            #print e
            thetarfile = tfile + e + filename
            page = urlopen(thetarfile)
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            #print gzip_f
            vcf = parser.Reader(gzip_f, "r")
            vcf._stream = thetarfile
            resulting_vcf.append(vcf)
            print (resulting_vcf)
        return resulting_vcf



def thousandgenomes_geo(**kwargs):
    '''This method enables to select VCF files corresponding to Arabidopsis Thaliana strains living at chosen
    longitude and / or latitude.

    Longitude and latitude should be provided as an interval (int1, int2).

    Method returns list of VCF.Reader objects created from selected file streams. '''

    longitude = kwargs.get('longitude', None)
    latitude = kwargs.get('latitude', None)
    f = open("thaliana_strains.csv", "r+")
    f.readline()
    e_list = []
    if longitude and not latitude:
        for line in f:
            l = line.split(',')
            if l[5] != "":
                if longitude[0] <= float(l[5]) <= longitude[1]:
                    e_list.append(l[0])

    elif latitude and not longitude:
        for line in f:
            l = line.split(',')
            if l[4]!="":
                if latitude[0] <= float(l[4]) <= latitude[1]:
                    e_list.append(l[0])

    elif latitude and longitude:
        for line in f:
            l = line.split(',')
            if l[4]!= "" and l[5] != "":
                if (latitude[0] <= float(l[4]) <= latitude[1] and longitude[0] <= float(l[5]) <= longitude[1]):
                    e_list.append(l[0])
    tfile = "http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/"
    filename = "_snp_short_indel_with_quality_reference.vcf.gz"
    if sys.version > '3':
        resulting_vcf = []
        for e in e_list:
            thetarfile = tfile + e + filename
            page = urlopen(thetarfile)
            gzip_f = gzip.GzipFile(mode='rb', fileobj=page)
            reader = codecs.getreader("utf-8")
            contents = reader(gzip_f)
            vcf = parser.Reader(contents, 'r')
            vcf._stream = thetarfile
            resulting_vcf.append(vcf)
            print(resulting_vcf)
        return resulting_vcf
    elif sys.version < '3':
        resulting_vcf = []
        for e in e_list:
            # print e
            thetarfile = tfile + e + filename
            page = urlopen(thetarfile)
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            # print gzip_f
            vcf = parser.Reader(gzip_f, "r")
            vcf._stream = thetarfile
            resulting_vcf.append(vcf)
            print(resulting_vcf)
        return resulting_vcf

def download(vcf_reader, path_filename):
    '''Downloads VCF file corresponding to stream on which provided VCF.Reader object is initialized.

    VCF.Reader object as well as download directory and filename are required.'''

    urlretrieve(vcf_reader._stream, path_filename)
    return









