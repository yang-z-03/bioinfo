
# utilities of reading fastq files from raw or gzipped archives.

from ansi import error

# read only the sequence partition of the file, and assign it with the
# corresponding name. ignoring the quality parts.

def readfq_seqs(fname, encoding = 'utf-8'):
    
    with open(fname, 'r', encoding = encoding) as f:
        return readfq_seqs_fp(f)
    
    pass

def readfq_seqs_gzipped(fname, encoding = 'utf-8'):

    import gzip
    import io

    with gzip.open(fname, 'rb') as f:
        with io.TextIOWrapper(f, encoding = encoding) as enc:
            return readfq_seqs_fp(enc)

    pass

def readfq_seqs_fp(fp):

    line = fp.readline().replace('\n', '')
    lineno = 0

    name = ''
    lines = {
        'names': [],
        'seqs': []
    }

    while line:

        # process the lines

        if lineno % 4 == 0:
            lines['names'] += [line]
        
        elif lineno % 4 == 1:
            lines['seqs'] += [line]

        line = fp.readline().replace('\n', '')
        lineno += 1
    
    if len(lines['names']) != len(lines['seqs']):
        error('inconsistant length of name and seqs in fastq file. ' + 
              'possibly unexpected end of file.')
    
    return lines