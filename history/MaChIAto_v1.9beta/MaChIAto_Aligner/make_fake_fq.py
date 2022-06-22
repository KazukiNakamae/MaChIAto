import sys
import re

fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
re_fa_header = re.compile('^>(.+)')

fa_header = ''
fa_seq    = ''
qa_header = '+\n'
qa_seq    = ''

if 2 > sum([1 for _ in open(fasta_file, 'r')]):
    with open(fastq_file, 'a') as fh2:
        fh2.write('\n')
    sys.exit()

with open(fasta_file, 'r') as fh:
    # import pdb; pdb.set_trace()
    for buf in fh:
        buf.rstrip('rn')
        m = re_fa_header.match(buf)
        if m:
            fa_header = '@' + m.group(1) + '\n'
            fa_seq    = ''
            qa_seq    = ''
        else:
            fa_seq = fa_seq + buf
            qa_seq = qa_seq + ''.join(['I'] * (len(buf) - 1)) + '\n'
            if fa_header != '':
                # print(fa_header, fa_seq, qa_header, qa_seq)
                with open(fastq_file, 'a') as fh2:
                    fh2.write(fa_header + fa_seq + qa_header + qa_seq)
sys.exit()

# print(fa_header + fa_seq + qa_header + qa_seq, )