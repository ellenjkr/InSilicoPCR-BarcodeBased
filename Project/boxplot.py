from Bio import SeqIO
import matplotlib.pyplot as plt
import os

# Lista de diretórios contendo os arquivos FASTA
input_directories = ["ANI1", "ANI2", "ANI3", "BAC", "VEG1", "VEG2", "COM"]

# Lista para armazenar as contagens de reads para cada pasta
read_counts = []
max_reads = 0
# Iterar pelos diretórios
for input_directory in input_directories:
    folder_read_counts = []

    # Iterar pelos arquivos no diretório
    for filename in os.listdir(input_directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(input_directory, filename)
            with open(filepath, "r") as fasta_file:
                records = list(SeqIO.parse(fasta_file, "fasta"))
                if len(records) > max_reads:
                    max_reads = len(records)
                folder_read_counts.append(len(records))
    if input_directory == "VEG2":
        folder_read_counts.append(0)

    read_counts.append(folder_read_counts)

# Crie um gráfico boxplot para cada pasta
plt.figure(figsize=(18, 6))
plt.boxplot(read_counts, labels=["ANI1", "ANI2", "ANI3", "BAC", "VEG1", "VEG2", "COM"], widths=0.5)  # Ajuste a largura aqui
plt.title("Distribuição de leituras  por conjunto")
plt.ylabel("Número de leituras")
plt.xlabel("Conjunto")
plt.yticks(range(0, max_reads + 1, 200))

plt.yscale('log')
import matplotlib.ticker as ticker
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.gca().set_yticks([5, 10, 20, 50, 100, 200, 500, 800, 2000, 3000, 5000, 8000, 10000])
plt.grid(axis='y')
plt.gca().set_axisbelow(True)
plt.show()
