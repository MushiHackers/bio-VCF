from Bio import VCF
reader = VCF.PhasedReader('Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_yri.D.phased')
##to Ci daje caly plik wczytany jako obiekt reader
#print(reader.filedata)
##to informacje o pliku jesli sa dostepne
#print(reader.filedata['region'])
##to slownik wiec tak mozna sie odwolywac
#print(reader.haplotypes)
##to jest lista haplotypow
#print(reader.haplotypes[0].name,reader.haplotypes[0].is_transmitted)
##wezmy sobie obejrzyjmy ten pierwszy
#record = reader.next()
##bierzemy pierwszy rekord
#print(record)
#print(record.rsID)
#print(record.pos)
#pytanie czy mam mu dodac opcje zeby bral tez chromosomy tu wypisywal?
#print(record.samples)
##to jest lista probek (zasady dla danego snp)
'''for s in record.samples:
    print(str(s))'''
##tak sa wypisywanie probki
#print(record.samples[0].exists,record.samples[0].is_unresolved,record.samples[0].nucleotide,record.samples[0].haplotype)
## wszystkie informacje dostepne dla probki
##
#reader2 = VCF.PhasedReader('Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_asw.unr.phased.gz')
##radzi sobie tez z gzipowanymi
##streamow jeszcze nie robilam i nie testowalam
#reader.get_snp_with_specific_id('rs10400036')
##sprawdza, czy SNP o danym rsID jest w pliku - jezeli jest, to zwraca caly record z nim, jezeli nie - informacje, ze nie ma takiego w pliku
plik = open('plikphased.phased','w')
writer = VCF.PhasedWriter(plik, reader)
for record in reader.fetch(region='191761-112976029'):
    writer.write_record(record)
#writer.flush()
writer.close()

plikvcf = open('Tests/VCF/chr13.vcf')

r = VCF.VCFReader(plikvcf)

for x in r.fetch(13,())
