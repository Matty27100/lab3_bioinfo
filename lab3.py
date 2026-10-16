import matplotlib.pyplot as plt

fasta_path = "GCF_000843065.1_ViralProj14272_genomic.fna"  # Downloaded Corndog genome.

# Standard RNA genetic code:
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
genetic_code = {
    "UUU": "Phenylalanine", "UUC": "Phenylalanine", "UUA": "Leucine", "UUG": "Leucine",
    "UCU": "Serine", "UCC": "Serine", "UCA": "Serine", "UCG": "Serine",
    "UAU": "Tyrosine", "UAC": "Tyrosine", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cysteine", "UGC": "Cysteine", "UGA": "Stop", "UGG": "Tryptophan",
    "CUU": "Leucine", "CUC": "Leucine", "CUA": "Leucine", "CUG": "Leucine",
    "CCU": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "CAU": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "CGU": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AUU": "Isoleucine", "AUC": "Isoleucine", "AUA": "Isoleucine", "AUG": "Methionine",
    "ACU": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "AAU": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "AGU": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GUU": "Valine", "GUC": "Valine", "GUA": "Valine", "GUG": "Valine",
    "GCU": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "GAU": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGU": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}

# Read the sequence, ignoring the FASTA header and joining sequence lines.
sequence_parts = []
header_count = 0
try:
    with open(fasta_path, encoding="utf-8-sig") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                header_count += 1
                if header_count > 1:
                    raise SystemExit("Please use a FASTA file containing one sequence.")
            else:
                sequence_parts.append(line.upper())
except OSError as error:
    raise SystemExit(f"Could not read the FASTA file: {error}")

sequence = "".join(sequence_parts)

# Slide a window of 3 bases forward by 1 base each time.
frequencies = {}
for i in range(len(sequence) - 2):
    codon = sequence[i:i + 3]
    if all(base in "ACGT" for base in codon):
        frequencies[codon] = frequencies.get(codon, 0) + 1

if not frequencies:
    raise SystemExit("No valid DNA triplets were found.")

top_10 = sorted(frequencies.items(), key=lambda item: item[1], reverse=True)[:10]
total = sum(frequencies.values())

print(f"\nSequence length: {len(sequence)} bases")
print(f"Valid triplets counted: {total}")
print("\nTop 10 most frequent codons:")
print(f"{'DNA':<6} {'RNA':<6} {'Count':>10} {'Frequency':>10}  Amino acid")
for codon, count in top_10:
    rna_codon = codon.replace("T", "U")
    amino_acid = genetic_code[rna_codon]
    frequency = count / total * 100
    print(f"{codon:<6} {rna_codon:<6} {count:>10} {frequency:>9.2f}%  {amino_acid}")

# Display the top 10 codons in a bar chart.
codons = [codon for codon, count in top_10]
counts = [count for codon, count in top_10]

plt.figure(figsize=(10, 5))
bars = plt.bar(codons, counts, color="steelblue")
plt.bar_label(bars, padding=3)
plt.title("Top 10 most frequent codons")
plt.xlabel("DNA codon")
plt.ylabel("Count")
plt.ylim(0, max(counts) * 1.15)
plt.tight_layout()
plt.show()
