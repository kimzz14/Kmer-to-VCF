#version Kmer-to-VCF 1.0.0
################################################################################
from optparse import OptionParser
import sys

#option parser
parser = OptionParser(usage="""Run kmer-to-vcf.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-o","--output",action = 'store',type = 'string',dest = 'OUTPUT',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.OUTPUT == None:
    print('Basic usage')
    print('')
    print('     python kmer-to-vcf.py -i kmc_files -o result.vcf')
    print('')
    print('Common options')
    print('')
    print('     -i      list of kmc files')
    print('     -o      output file name')
    sys.exit()

class KMC_READER:
    def __init__(self, fileName, readN):
        self.fin = open(fileName)
        self.hasNext = True
        self.readN = readN
        self.kmerIDX = self.readN
        self.next_line()

    def next_line(self):
        line = self.fin.readline()
        if not line:
            self.hasNext = False
            return
        self.cSequence, self.cCount = line.rstrip('\n').split('\t')
        self.cKmerIDX = self.quaternary(self.cSequence)
    def has_next(self):
        return self.hasNext
    def next(self):
        data_DICT = {}
        while self.has_next():
            if self.cKmerIDX > self.kmerIDX:
                self.kmerIDX += self.readN
                break
            #data_DICT[self.cSequence] = (self.cCount, self.cKmerIDX)
            data_DICT[self.cSequence] = self.cCount
            self.next_line()
        return data_DICT

    def quaternary(self, sequence):
        nucl_DICT = {'A':0, 'C': 1, 'T': 2, 'G': 3}
        count = 0
        for idx, nucl in enumerate(sequence[::-1]):
            count += nucl_DICT[nucl] * pow(4, idx)
        return count
################################################################################
sample_LIST = []
sample_DICT = {}

fin = open(opt.INPUT)
for line in fin:
    sample = line.rstrip('\n')
    sample_LIST += [sample]
fin.close()

for sample in sample_LIST:
    sample_DICT[sample] = KMC_READER(sample, 100000000)
################################################################################
fout = open(opt.OUTPUT, 'w')
fout.write("##fileformat='VCFv4" + '\n')
fout.write('##FORMAT=<ID=GT, Number=1, Type=String, Description="Genotype"' + '\n')
fout.write('##contig=<ID=KMER, Description="Virtual chromosome for VCF"' + '\n')
fout.write('#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_LIST) + '\n')
pos = 0
while True:
    tmp_DICT = {}
    isAllEnd = False
    for sample in sample_LIST:
        if sample_DICT[sample].has_next() == True:
            isAllEnd = True
            break
    if isAllEnd == False: break
    key_LIST = []
    for sample in sample_LIST:
        tmp_DICT[sample] = sample_DICT[sample].next()
        key_LIST += tmp_DICT[sample].keys()
    key_LIST = sorted(set(key_LIST))

    for key in key_LIST:
        genotype_LIST = []
        for sample in sample_LIST:
            if key in tmp_DICT[sample]:
                genotype_LIST += ['0/0']
            else:
                genotype_LIST += ['1/1']
        if len(set(genotype_LIST)) == 1:
            continue
        pos += 1
        context = ['KMER', str(pos), key, 'A', 'T', '.', '.', '.', 'GT']
        fout.write('\t'.join(context + genotype_LIST) + '\n')
fout.close()
################################################################################
print('Done')
