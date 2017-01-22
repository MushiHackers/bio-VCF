from Bio import VCF
from Bio.VCF.parser import _Filter, Reader
import sys


class Filters:

    def __init__(self, arguments):
        self.arguments=arguments

    def CHR(self, record,val=None): #Chromosom number(s)
        try:
            for i in range(len(val)):
                val[i]=str(val[i])
            if record.CHROM.split('chr')[-1] in val:
                return True
            else:
                return False
        except:
            sys.exit("check if you gave list of chromosomes names")

    def POS(self, record,val=None): #Chromosom number(s)
        try:
            if record.POS in val:
                return True
            else:
                return False
        except:
            sys.exit("check if you gave list of position number")

    def SNP(self, record,val=None): #SnpOnly
        if record.is_snp==True:
            return True
        else:
            return False

    def AVG(self,record,val): #AvgDepthPerSample
        try:
            avgcov = float(record.INFO['DP']) / len(record.samples)
            if avgcov < val:
                return True #avgcov
            else:
                return False
        except:
            sys.exit("can't use filter AVG, no information about DP")

    def DPS(self,record,val=None): #DepthPerSample
        try:
            if record.is_indel==True:
                return False

            mindepth = min([sam['DP'] for sam in record.samples])
            if mindepth < val:
                return True #mindepth
            else:
                return False
        except:
            sys.exit("can't use filter DPS, no information about DP")

    def VGQ(self,record,val): #VariantGenotypeQuality
        try:
            if not record.is_monomorphic:
                vgq = max([x['GQ'] for x in record if x.is_variant])
                if vgq < val:
                    return True #vgq
                else:
                    return False
            else:
                return False
        except:
            sys.exit("can't use filter VGQ, no information about GQ")

    def SQ(self,record,val): #SiteQuality
        if record.QUAL < val:
            return True #record.QUAL
        else:
            return False

    def recognition(self, name, record,argument):
        if name=="SNP" or name=='snp' or name=='Snp':
            return self.SNP(record,argument)

        if name=="AVG" or name=='avg' or name=='Avg':
            return self.AVG(record,argument)

        if name=="DPS" or name=='dps' or name=='Dps':
            return self.DPS(record,argument)

        if name=="VGQ" or name=='vgq' or name=='Vgq':
            return self.VGQ(record,argument)

        if name=="SQ" or name=='sq' or name=='Sq':
            return self.SQ(record,argument)

        if name=="CHR" or name=='chr' or name=='Chr':
            return self.CHR(record,argument)

        if name=="POS" or name=='pos' or name=='Pos':
            return self.POS(record,argument)
        else:
            sys.exit("unknown filter name %s" % (name))


    def __call__(self,record):
        for i in self.arguments:
            if self.recognition(i, record, self.arguments[i])!=True:
                return False
        return True

class Base:

    def __init__(self, vcf_file, output_file=None):
        self.reader=Reader(open(vcf_file, 'r'))
        self.output_file=output_file

    def filtering(self, filter_dict):

        if self.output_file!=None:
            f=open(self.output_file,'w')
        filt=Filters(filter_dict)
        for record in self.reader:
            if filt(record)==True:
                if self.output_file!=None:
                    f.write(record+'\n')
                print(record)
        if self.output_file!=None:
            f.close()
        return True


