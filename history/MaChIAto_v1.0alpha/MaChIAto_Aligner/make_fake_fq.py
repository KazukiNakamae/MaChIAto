import sys
import re

fasta_file = sys.argv[1]
re_fa_header = re.compile('^>(.+)')

fa_header = ''
fa_seq    = ''
qa_header = '+\n'
qa_seq    = ''

with open(fasta_file, 'r') as fh:
    for buf in fh:
        buf.rstrip('rn')
        m = re_fa_header.match(buf)
        if m:
            if fa_header != '':
                print(fa_header, fa_seq, qa_header, qa_seq)
            fa_header = '@' + m.group(1) + '\n'
            fa_seq    = ''
            qa_seq    = ''
        else:
            fa_seq = fa_seq + buf
            qa_seq = qa_seq + ''.join(['I'] * (len(buf) - 1))

print(fa_header, fa_seq, qa_header, qa_seq)