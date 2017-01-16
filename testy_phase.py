import sys
import re
#import pybedtools
from Bio import VCF
reader = VCF.PhasedReader('Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_yri.D.phased')
reader2 = VCF.PhasedReader('Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_yri.D.phased')

'''rec1 = reader.next()
rec2 = reader2.next()'''

'''print(rec1==rec2)
sys.exit()'''
##to Ci daje caly plik wczytany jako obiekt reader
# print(reader.filedata)
# ##to informacje o pliku jesli sa dostepne
# print(reader.filedata['region'])
# ##to slownik wiec tak mozna sie odwolywac
# print(reader.haplotypes)
# ##to jest lista haplotypow
# print(reader.haplotypes[0].name,reader.haplotypes[0].is_transmitted)
# ##wezmy sobie obejrzyjmy ten pierwszy
# record = reader.next()
# ##bierzemy pierwszy rekord
# print(record)
# print(record.rsID)
# print(record.pos)
# #pytanie czy mam mu dodac opcje zeby bral tez chromosomy tu wypisywal?
# print(record.samples)
# ##to jest lista probek (zasady dla danego snp)
# for s in record.samples:
#     print(str(s))
# ##tak sa wypisywanie probki
# print(record.samples[0].exists,record.samples[0].is_unresolved,record.samples[0].nucleotide,record.samples[0].haplotype)
## wszystkie informacje dostepne dla probki
##
#reader2 = VCF.PhasedReader('Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_asw.unr.phased.gz')
##radzi sobie tez z gzipowanymi
##streamow jeszcze nie robilam i nie testowalam
#print reader.get_snp_with_specific_id('rs2066314')
##sprawdza, czy SNP o danym rsID jest w pliku - jezeli jest, to zwraca caly record z nim, jezeli nie - informacje, ze nie ma takiego w pliku
#print reader.get_snp_within_range(418076, 504032)
##wypisuje wszystkie snpy w danym zakresie (zakres: >= i <)
#print reader.get_snp_with_specific_sample('NA18855_B', 'T')
##wypisuje wszystkie snpy z dana probka
#print reader.get_samples_in_specific_hap('NA18855_NA18856_A')
##wypisuje wszystkie snpy w danym haplotypie
#reader.get_sample_from_specific_snp('rs11252546', 'NA18855_NA18856_A')
##wypisuje zasade w danym snpie w danym haplotypie
#plik = open('plikphased.phased','w')
#writer = VCF.PhasedWriter(plik, reader)
#for record in reader.fetch(region='191761-112976029'):
#    writer.write_record(record)
#writer.flush()
#writer.close()



'''print('-a-a-a-a-a-a-')
plik = open('pliczko.phased','w')
plikvcf = open('Tests/VCF/chr10.vcf')
writer = VCF.PhasedWriter(plik, reader)'''
# read = VCF.Reader(plikvcf,prepend_chr=True)
# read = read.fetch(chrom='1',verbose=False, vcf='temp')
#
# for v in read:
#
#     if isinstance(read, pybedtools.bedtool.BedTool):
#         if len(v[3]) > 1:
#             continue
#         else:
#             _alt_pattern = re.compile('[\[\]]')
#             alter = v[4].split(',')
#             for alt in alter:
#                 if _alt_pattern.search(alt) is not None:
#                     continue
#                 elif alt[0] == '.' and len(alt) > 1:
#                     continue
#                 elif alt[-1] == '.' and len(alt) > 1:
#                     continue
#                 elif alt[0] == "<" and alt[-1] == ">":
#                     continue
#                 elif alt not in ['A', 'C', 'G', 'T', 'N', '*']:
#                     continue
#                 else:
#                     print('je!')
#                     print(v)
#     else:
#         if v.is_snp:
#             print('je2!')
#             print(v)
#
# sys.exit()


'''newr = reader.fetch(fsock=plikvcf, verbose=False)
for record in newr:
    writer.write_record(record)'''

#newr = reader.fetch(region='191761-112976029')
#print(newr.filedata)
##to informacje o pliku jesli sa dostepne
#print(newr.filedata['region'])
##to slownik wiec tak mozna sie odwolywac
#print(newr.haplotypes)
##to jest lista haplotypow
#print(newr.haplotypes[0].name,newr.haplotypes[0].is_transmitted)
##wezmy sobie obejrzyjmy ten pierwszy
#record = newr.next()
#record = newr.next()
##bierzemy pierwszy rekord
#print(record)
#print(record.rsID)
#print(record.pos)
#pytanie czy mam mu dodac opcje zeby bral tez chromosomy tu wypisywal?
#print(record.samples)
##to jest lista probek (zasady dla danego snp)
#for s in record.samples:
#    print(str(s))
##tak sa wypisywanie probki
#print(record.samples[1].exists,record.samples[1].is_unresolved,record.samples[1].nucleotide,record.samples[1].haplotype,record.samples[1].is_not_matching_snp)
## wszystkie informacje dostepne dla probki
sys.exit()

'''print('--------------')

print(reader.filedata)
##to informacje o pliku jesli sa dostepne
print(reader.filedata['region'])
##to slownik wiec tak mozna sie odwolywac
print(reader.haplotypes)
##to jest lista haplotypow
print(reader.haplotypes[0].name,reader.haplotypes[0].is_transmitted)
##wezmy sobie obejrzyjmy ten pierwszy
record = reader.next()
##bierzemy pierwszy rekord
print(record)
print(record.rsID)
print(record.pos)
#pytanie czy mam mu dodac opcje zeby bral tez chromosomy tu wypisywal?
print(record.samples)
##to jest lista probek (zasady dla danego snp)
for s in record.samples:
    print(str(s))
##tak sa wypisywanie probki
print(record.samples[0].exists,record.samples[0].is_unresolved,record.samples[0].nucleotide,record.samples[0].haplotype)'''
