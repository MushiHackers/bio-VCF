import gzip, io, codecs, sys
from Bio.VCF import parser

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
    from urllib.request import urlretrieve
except ImportError:
    # Fall back to Python 2's urllib2 and urllib
    from urllib2 import urlopen
    from urllib import urlretrieve


def thousandgenomes(name = None, ecotype = None, ecnumber = None, country = None, longitude = None, latitude = None):
    '''This method enables to search through '1001 Genomes' database for VCF file corresponding to selected
        Arabidopsis Thaliana strain, origin country or longitude and/or latitude where Arabidopsis Thaliana live.

        Country name must be chosen from: "USA", "FRA", "CZE", "AUT", "KGZ", "TJK", "SWE", "UK", "GER", "KAZ",
        "BEL", "CPV", "ESP", "RUS", "NED", "FIN", "SUI", "ITA", "IRL", "POR", "EST", "DEN", "IND", "LTU", "JPN", "POL", "NOR",
        "CAN", "UKR", "AZE", "GEO", "ARM", "MAR", "CRO", "BUL", "GRC", "SVK", "ROU", "UZB", "SRB", "CHN", "IRN", "LBN", "MAR",
        "AFG".

        Longitude and latitude should be provided as an interval - longitude = (int1, int2), latitude = (int1,int2).

        Search is made based on provided strain name / Eco type / EC number / country name / longitude / latitude.

        Method returns VCF.Reader object of selected VCF file stream from 1001 Genomes Database.'''

    if name:
        return _thousandgenomes_more(name = name)
    if ecotype:
        return _thousandgenomes_more(ecotype = ecotype)
    if ecnumber:
        return _thousandgenomes_more(ecnumber=ecnumber)
    if country:
        return _thousandgenomes_more(country=country)
    if longitude and not latitude:
        return _thousandgenomes_geo(longitude=longitude)
    if latitude and not longitude:
        return _thousandgenomes_geo(latitude=latitude)
    if latitude and longitude:
        return _thousandgenomes_geo(latitude=latitude, longitude=longitude)
    else:
        raise Exception('You must provide at least one argument')


def _thousandgenomes_more(file='thaliana_strains.csv', name = None, ecotype = None, ecnumber = None, country = None):
    '''Function used to search through '1001 Genomes Dabatase' in thousandgenomes() method.'''

    if not (ecotype or name or ecnumber or country):
        raise Exception('You must provide strain name / eco type / EC number or country')
    if not ecotype:
        if name:
            pattern = ','+str(name)+','
            f = open(file,"r+")
            for line in f:
                if pattern in line:
                    l=line.split(',')
                    ecotype = l[0]
            if not ecotype:
                raise Exception('Provided name does not refer to any strain name in our database')
        if ecnumber:
            pattern = ','+str(ecnumber)+','
            f = open(file,"r+")
            for line in f:
                if pattern in line:
                    l = line.split(',')
                    ecotype = l[0]
            if not ecotype:
                raise Exception('Provided EC number does not refer to any strain in out database')
        if country:
            f = open(file, "r+")
            pattern = ',' + str(country) + ','
            e_list = []
            for line in f:
                if pattern in line:
                    l = line.split(',')
                    e_list.append(l[0])

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
        else: # sys.version < '3':
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            vcf = parser.Reader(gzip_f,"r")
            vcf._stream = thetarfile
            return vcf
    if e_list:
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
                #print(resulting_vcf)
            return resulting_vcf
        else: # sys.version < '3':
            resulting_vcf = []
            for e in e_list:
                # print e
                thetarfile = tfile + e + filename
                page = urlopen(thetarfile)
                gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
                vcf = parser.Reader(gzip_f, "r")
                vcf._stream = thetarfile
                resulting_vcf.append(vcf)
                print(resulting_vcf)
            return resulting_vcf


def _thousandgenomes_geo(file='thaliana_strains.csv', longitude = None, latitude = None):
    '''Function used to search through '1001 Genomes Dabatase' in thousandgenomes() method.'''

    f = open(file, "r+")
    f.readline()
    e_list = []
    if (longitude and not latitude):
        for line in f:
            l = line.split(',')
            if l[5] != "":
                if longitude[0] <= float(l[5]) <= longitude[1]:
                    e_list.append(l[0])

    elif (latitude and not longitude):
        for line in f:
            l = line.split(',')
            if l[4]!="":
                if latitude[0] <= float(l[4]) <= latitude[1]:
                    e_list.append(l[0])

    elif (latitude and longitude):
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
            #print(resulting_vcf)
        return resulting_vcf
    else: # sys.version < '3':
        resulting_vcf = []
        print (resulting_vcf)
        for e in e_list:
            # print e
            thetarfile = tfile + e + filename
            page = urlopen(thetarfile)
            gzip_f = gzip.GzipFile(fileobj=io.BytesIO(page.read()))
            # print gzip_f
            vcf = parser.Reader(gzip_f, "r")
            vcf._stream = thetarfile
            resulting_vcf.append(vcf)
            #print(resulting_vcf)
        return resulting_vcf

def download(vcf_reader, path_filename):
    '''Downloads VCF file corresponding to stream on which provided VCF.Reader object is initialized.

    VCF.Reader object as well as download directory and filename are required.'''

    urlretrieve(vcf_reader._stream, path_filename)
    return









