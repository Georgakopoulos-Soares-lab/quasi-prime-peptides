import itertools
import time
import json
import os
import logging


class PeptideFrequencies:
    """
    Calculate the frequency of peptides in a category (e.g. in bacteria)

    Input a list of files containing the peptides in each species proteome separated by newlines. Calculates the frequency of each peptide in the category
    from 0 to 1. Does not check the uniqueness of peptides in the input files - if an input file contains duplicates of the same peptide frequencies higher than
    1 can be returned.

    An example run can be found at the bottom of the file.
    """
    def __init__(self, amino_acids_list: list[str], peptide_length: int, output_folder: str):
        self.amino_acids_list = amino_acids_list
        self.peptide_length = peptide_length
        self.output_folder = output_folder
        self.timestamp = time.strftime('%Y%m%d_%H%M-', time.localtime())
        if not os.path.exists(self.output_folder+"/logs"):
            os.makedirs(self.output_folder+"logs")
        logging.basicConfig(filename=self.output_folder+"logs/"+self.timestamp+"peptide_frequencies.log", encoding='utf-8', level=logging.DEBUG)

    def analyze(self, files):
        self.count_peptides(files)
        self.calculate_frequencies(len(files))

    def count_peptides(self, files: list[str]):
        # Generate a dictionary with all possible kmers from alphabet AMINO_ACIDS LIST and length PEPTIDE_LENGTH
        self.peptides_dict: dict[str, float] = dict.fromkeys((''.join(p) for p in itertools.product(self.amino_acids_list, repeat=self.peptide_length)), 0)

        logging.info(f"Analyzing peptide frequencies in the following files: {files}")

        for file in files:
            with open(file, "r") as f:
                for line in f:
                    try:
                        self.peptides_dict[line.strip()] += 1
                    except KeyError:
                        logging.warning(f"Peptide {line.strip()} was not found in the dictionary")

    def calculate_frequencies(self, num_files: int):
        for peptide in self.peptides_dict:
            self.peptides_dict[peptide] = self.peptides_dict[peptide] / num_files

    def save_as_json(self, prefix: str):
        with open(self.output_folder + self.timestamp + prefix + ".json", "w") as f:
            json.dump(self.peptides_dict, f)


if __name__ == "__main__":
    import glob

    AMINO_ACIDS_LIST = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    PEPTIDE_LENGTH = 5
    OUTPUT_FOLDER = "output/peptide_frequencies/"
    PF = PeptideFrequencies(AMINO_ACIDS_LIST, PEPTIDE_LENGTH, OUTPUT_FOLDER)

    categories = ["bacteria", "archaea", "viruses", "eukaryota"]
    for category in categories:
        files = glob.glob("test_files/*" + category + "*5mers")
        PF.analyze(files)
        PF.save_as_json(category)


