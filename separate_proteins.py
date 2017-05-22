with open('test.faa', 'r') as fin, open('test.out.faa', 'w') as fout:
    for line in fin:
        if '>' in line:
            hdr = line
        else:
            protein_seqs = line.split('*')
            for s in protein_seqs[:-1]:
                fout.write(hdr)
                fout.write(s + '\n')
            fout.write(hdr)
            fout.write(s)