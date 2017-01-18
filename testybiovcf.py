from __future__ import print_function
from Bio import VCF

print(1)
vcf_reader = VCF.Reader(open('/home/hania/Desktop/apb/bio-VCF/Tests/VCF/example-4.0.vcf', 'r'))
for record in vcf_reader:
     print(record)

print(2)
vcf_reader = VCF.Reader(open('/home/hania/Desktop/apb/bio-VCF/Tests/VCF/example-4.0.vcf', 'r'))
record = next(vcf_reader)
print(record.POS)
print(record.ALT)
print(record.INFO['AF'])

print(3)
print(record.num_called, record.call_rate, record.num_unknown)
print(record.num_hom_ref, record.num_het, record.num_hom_alt)
print(record.nucl_diversity, record.aaf, record.heterozygosity)
print(record.get_hets())
print(record.is_snp, record.is_indel, record.is_transition, record.is_deletion)
print(record.var_type, record.var_subtype)
print(record.is_monomorphic)

print(4)
record = next(vcf_reader)
for sample in record.samples:
    print(sample['GT'])
print(record.genotype('NA00001')['GT'])

print(5)
call = record.genotype('NA00001')
print(call.site)
print(call.sample)
print(call.data)

print(6)
print(call.called, call.gt_type, call.gt_bases, call.phased)

print(7)
print(vcf_reader.metadata['fileDate'])
print(vcf_reader.samples)
print(vcf_reader.filters)
print(vcf_reader.infos['AA'].desc)

print(8)
reader = VCF.Reader(open('/home/hania/Desktop/apb/bio-VCF/Tests/VCF/example-4.1-bnd.vcf'))
_ = next(reader); row = next(reader)
print(row)
bnd = row.ALT[0]
print(bnd.withinMainAssembly, bnd.orientation, bnd.remoteOrientation, bnd.connectingSequence)

print(9)
print("deprecated-fetch")
"""to nie jest potrzebne - bo to fetch; do zmiany"""
"""vcf_reader = VCF.Reader(filename='/home/hania/Desktop/apb/bio-VCF/Tests/VCF/tb.vcf.gz')
for record in vcf_reader.fetch('20', 1110695, 1230237):
    print(record)
print(vcf_reader.fetch('4', 10, 20))"""

print(10)
vcf_reader = VCF.Reader(filename='/home/hania/Desktop/apb/bio-VCF/Tests/VCF/tb.vcf.gz')
vcf_writer = VCF.Writer(open('/home/hania/Desktop/apb/bio-VCF/test10', 'w'), vcf_reader)
for record in vcf_reader:
    vcf_writer.write_record(record)

print('FILTERS')
"""do przemyslenia"""
print('todo')

print(11)
print(VCF.trim_common_suffix('TATATATA', 'TATATA'))
print(VCF.trim_common_suffix('ACCCCC', 'ACCCCCCCC', 'ACCCCCCC', 'ACCCCCCCCC'))

t = VCF.PhasedReader(filename='./Tests/VCF/hapmap3_r2_b36_fwd.consensus.qc.poly.chr10_yri.D.phased.gz')
t2 = t.fetch(fsock=open('./Tests/VCF/chr10.vcf','r'),vcf='testowy')
