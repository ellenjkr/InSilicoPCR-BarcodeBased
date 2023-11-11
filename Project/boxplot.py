from Bio import SeqIO
import matplotlib.pyplot as plt
import os


def run(config, primer_names):
    # Lista de diretórios contendo os arquivos FASTA


    input_directories = [item for item in os.listdir(config['OUTPUT_PATH']) if os.path.isdir(config['OUTPUT_PATH'])] 

    # Lista para armazenar as contagens de reads para cada pasta
    read_counts = []
    max_reads = 0
    os.chdir(config['OUTPUT_PATH'])
    # Iterar pelos diretórios
    for input_directory in input_directories:
        folder_read_counts = []

        # Iterar pelos arquivos no diretório
        for filename in os.listdir(input_directory):
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                file_basename = filename.replace('.fasta', '')
                file_basename = file_basename.replace('.fa', '')

                if file_basename not in primer_names and config['HAS_BARCODES'] == 'true' or config['HAS_BARCODES'] == 'false':
                    filepath = os.path.join(input_directory, filename)
                    with open(filepath, "r") as fasta_file:
                        records = list(SeqIO.parse(fasta_file, "fasta"))
                        if len(records) > max_reads:
                            max_reads = len(records)
                        folder_read_counts.append(len(records))


        read_counts.append(folder_read_counts)

    # Crie um gráfico boxplot para cada pasta
    plt.figure(figsize=(18, 6))
    plt.boxplot(read_counts, labels=input_directories, widths=0.5)  # Ajuste a largura aqui
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


if __name__ == '__main__':
    run({'OUTPUT_PATH': 'output/'}, 'true', ['BAC1'])