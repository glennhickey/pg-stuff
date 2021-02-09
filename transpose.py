import sys

def load(f):
    mat = []
    with open(f, 'r') as infile:
        for line in infile:
            toks = line.strip().split()
            if toks:
                mat.append(toks)
    return mat

for i, arg in enumerate(sys.argv[1:]):
    mat = load(arg)
    if i == 0:
        sys.stdout.write('name\t')
        header = [row[0] for row in mat]
        sys.stdout.write('\t'.join(header))
        sys.stdout.write('\n')
    name = arg.split('.')[0]
    vals = [name] + [row[1] for row in mat]
    sys.stdout.write('\t'.join(vals))
    sys.stdout.write('\n')
            
        
