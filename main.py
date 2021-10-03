from analysis import Sequence
import sys


def colored(seq):

    COLORS = {
        'A': '\033[92m' + 'A',
        'C': '\033[94m' + 'C',
        'G': '\033[93m' + 'G',
        'T': '\033[91m' + 'T',
        'U': '\033[91m' + 'U',
        'reset': '\033[0;0m'
    }

    result = ""

    for nuc in seq:
        result += COLORS[nuc] if nuc in COLORS else nuc
    return result + '\033[0;0m'


def read_file(fname):

    with open(fname, 'r') as f:
        return "".join([i.strip() for i in f.readlines() if not i.startswith('>')])



def main(filename=""):

    seq_data = read_file(filename) if filename else "" 


    sequence = Sequence(seq_data)


    terminal_template = f"""
    
    [!] Result Saved to output.txt


    Sequence: {colored(sequence.sequence[:60])} ...

    + Length: {len(sequence.sequence)}
    
    + Type: {sequence.type}

    + Nucleotide Frequency: {sequence.nc_frequency}
    
    + DNA/RNA Transcription: {colored(sequence.rna_transcription[:60]) + "..."}

    + DNA string + complement + reverse complement:
        

        5' {colored(sequence.sequence[:60])} 3'
           {'|' * 60}
        3' {colored(sequence.complement[:60])} 5' [Complement]
        5' {colored(sequence.rev_complement[:60])} 3' [Reverse Complement]


    + GC Content: {sequence.GC_content}%

    + Aminoacids Sequence: {sequence.aminoacid_translaion[:10]} ...

    + Possible protiens sequence:
    """

    print(terminal_template)
    for i in sequence.protiens:
        print(" ", i)
    print()

    sequence.write_output()



if __name__ == "__main__":

    if len(sys.argv) < 2:
        main()
    else:
        main(sys.argv[1])
