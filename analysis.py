import random
from collections import Counter



class Structure:

    NUCLEOTIDES = {
        "DNA": ['A', 'C', 'G', 'T'],
        "RNA": ['A', 'C', 'G', 'U']
    }

    DNA_Codons = {
        # 'M' - START, '_' - STOP
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "M",
        "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y",
        "TAA": "_", "TAG": "_", "TGA": "_"
    }

    RNA_Codons = {
        # 'M' - START, '_' - STOP
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UGU": "C", "UGC": "C",
        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "UUU": "F", "UUC": "F",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAU": "H", "CAC": "H",
        "AUA": "I", "AUU": "I", "AUC": "I",
        "AAA": "K", "AAG": "K",
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUG": "M",
        "AAU": "N", "AAC": "N",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UGG": "W",
        "UAU": "Y", "UAC": "Y",
        "UAA": "_", "UAG": "_", "UGA": "_"
    }

    ProtienDB = "https://blast.ncbi.nlm.nih.gov/blast/smartblast/smartBlast.cgi?CMD=submit&QUERY="

    def format_url(self, seq):
        """
        Search DNA/RNA databse for a sequence.
        """
        return self.ProtienDB + seq


class Sequence(Structure):
    """
    pass
    """

    def __init__(self, sequence="") -> None:

        self.sequence = sequence.upper() if sequence != "" else self.generate_sequence()

        self.type = self.__detect_type()

        self.nc_frequency = self.nucleotide_frequency()
        
        self.rna_transcription = self.transcription() if self.type == "DNA" else self.sequence
        
        self.complement = self.reverse_complement()[::-1]
        
        self.rev_complement = self.reverse_complement()
        
        self.GC_content = self.gc_content()

        self.aminoacid_translaion = self.translate()
        
        self.open_reading_frames = self.reading_frames()
        
        self.protiens = self.possible_proteins()

    def __detect_type(self):
        """
        Detect and validate a DNA/RNA sequence.
        """

        copy_sequence = set(self.sequence)

        if copy_sequence.issubset(self.NUCLEOTIDES["DNA"]):
            return "DNA"
        
        elif copy_sequence.issubset(self.NUCLEOTIDES["RNA"]):
            return "RNA"
        
        assert "Please provide a valid DNA/RNA sequence."
    
    def generate_sequence(self, length=60):
        """
        Generates a random DNA/RNA sequence.
        """

        sequence_type = random.choice(["DNA", "RNA"])
        return "".join([random.choice(self.NUCLEOTIDES[sequence_type]) for _ in range(length)]).upper()


    def nucleotide_frequency(self):
        """
        Coun nucleotides. 
        """
        
        return dict(Counter(self.sequence))
    
    def transcription(self):
        """
        DNA to RNA transcription.
        """
        
        # Replace thymine with uracil
        return self.sequence.replace('T', 'U')
    
    def reverse_complement(self):
        """
        Swapping adenine with thymine and guanine with cytosine.
        """

        mapping =  str.maketrans('ATCG', 'TAGC') if self.type == "DNA" else \
            str.maketrans('AUCG', 'UAGC')
        return self.sequence.translate(mapping)[::-1]

    def gc_content(self):
        """GC Content in a DNA/RNA sequence"""

        return round((self.sequence.count('C') + self.sequence.count('G')) / len(self.sequence) * 100)

    def translate(self, init_idx=0, seq=""):
        """
        Translates DNA/RNA sequence to amino acids sequence.
        """

        sequence = self.sequence if not seq else seq
        cordons = self.DNA_Codons if self.type == "DNA" else self.RNA_Codons
        return [cordons[sequence[i:i + 3]] for i in range(init_idx, len(sequence) - 2, 3)]


    def reading_frames(self):
        """Generate the six reading frames, including reverse complement"""
        
        result = [] 
        for _ in range(3):
            result.append(self.translate(_) )
        result.extend([
            self.translate(init_idx=_, seq=self.rev_complement) for _ in range(3)])
        return result
            

    def proteins_from_rf(self, aminoacid_seq):
        """
        Compute all possible proteins in an aminoacid seq.
        """
        
        current_prot, proteins = [], []
        
        for aa in aminoacid_seq:
            if aa == "_":
                if current_prot:
                    proteins.extend(current_prot)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    
    def possible_proteins(self):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""

        reading_frames = self.open_reading_frames
        possible_protiens = [
            self.proteins_from_rf(frame) for frame in reading_frames
        ]

        return sorted([self.format_url(protien) for _ in possible_protiens for protien in _], key=len, reverse=True)

    def write_output(self):
        """
        Write result to a file.
        """

        with open("output.txt", 'w') as f:

            file_template = f"""

            Sequence:
            
            {self.sequence}

            Length: {len(self.sequence)}
            
            Type: {self.type}
            
            GC Content: {self.GC_content}%

            Nucleotide Frequency:
            
            {self.nc_frequency}
            
            DNA/RNA Transcription:
            
            {self.rna_transcription}

            DNA complement:

            {self.complement}

            DNA Reverse Complement

            {self.rev_complement}

            Aminoacids Sequence:
            
            {self.aminoacid_translaion}

            Possible Protien Sequences:

            """
            f.write(file_template)
            for i in self.protiens:
                f.write('\t' + i + '\n')
