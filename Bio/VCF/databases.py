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

# ThousandGenomes
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

# dbSNP

def dbSNP_download(organism_taxon, chromosome = None):
    '''This method enable you to search VCF files corresponding to given species, taxon identifier and chromosome.
    
    You can find avaiable organisms and taxon identifiers in file organisms.txt
    
    If you want to check dbSNP for updates, or newly added organisms or VCF files use function check_VCF.
    
    Attention: human VCF file is not divided by chromosomes. You can download collective file for all chromosomes 
    but remember it's almost 350GB
    
    '''
    
    baseURL = "ftp://ftp.ncbi.nih.gov/snp/organisms/"
    if organism_taxon.split('_')[0] != 'human':
       vcf_url = baseURL + organism_taxon + '/' + 'VCF/' + 'vcf_chr_' + str(chromosome) + '.vcf.gz'
    else:
        search_url = baseURL + organism_taxon + '/' + 'VCF/'
        search = urlopen(search_url)
        for i in search:
            line = i.strip().split()[-1]
            name = line.split('_')
            if name[0] == 'All' and name[1].split('.')[-1] == 'gz' and name[1].split('.')[-2] == 'vcf':
                filename = line
        vcf_url = baseURL + organism_taxon + '/' + 'VCF/' + filename
     

    response = urlopen(vcf_url)
    
    if sys.version > '3':
        gzip_f = gzip.GzipFile(mode='rb', fileobj=response)
        reader = codecs.getreader("utf-8")
        contents = reader(gzip_f)
        vcf = parser.Reader(contents,"r")
        return vcf
    else: # sys.version < '3':
        gzip_f = gzip.GzipFile(fileobj=io.BytesIO(response.read()))
        vcf = parser.Reader(gzip_f,"r")
        return vcf

def check_VCF():
    
    '''This method updates the list of organisms in organisms.txt'''
    
    fullurl = "ftp://ftp.ncbi.nih.gov/snp/organisms/"
    database_list = urlopen(fullurl)
    organisms_list = []
    vcf_organism = []
    for i in database_list:
        a = i.strip().split()
        organisms_list.append(a[-1])
        
    for organism in organisms_list:
        if organism == 'FTP_MAINTENENCE_NOTICE.txt':
            pass
        else:
            organism_url = fullurl + organism + '/'
            check_organism = urlopen(organism_url)
            for line in check_organism:
                b = line.strip().split()
                if b[-1] == 'VCF':
                    #org = organism.split('_')
                    organism_vcf_url = organism_url + 'VCF/'
                    check_vcf = urlopen(organism_vcf_url)
                    file_list = []
                    for i in check_vcf:
                        file_list.append(i)
                    if len(file_list) > 2:
                        vcf_organism.append(organism)
                        
                    
    result = open('organisms.txt','w+')
    for record in vcf_organism:
        result.write(record + '\n')
    result.close()
 








