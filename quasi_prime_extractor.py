import os
import logging
import time


class Quasi_Prime_Extractor:
    """
    Use the extract function to identify all all kmers found in qp_files but not in background_files. Both variables need to be lists containing the paths of the
    files included in each category, where each file is a list of unique kmers present in that species separated by a newline (see also accompanying test files
    on github). If qp_files contains only one species, this results in identifying the quasi-primes of that species. If qp_files contains all the species of a
    taxonomy, this results in identifying the quasi-primes of that taxonomy.
    """
    def __init__(self, output_folder: str):
        self.output_folder = output_folder
        self.timestamp = time.strftime('%Y%m%d_%H%M', time.localtime())
        if not os.path.exists(self.output_folder + "/logs"):
            os.makedirs(self.output_folder + "logs")
        logging.basicConfig(filename=self.output_folder + "logs/" + self.timestamp + "peptide_frequencies.log", encoding='utf-8', level=logging.DEBUG)

    @staticmethod
    def read_kmers(file_list):
        """Generate a set containing the intersection of all kmer peptides in the input files
        Input files are assumed to be contain kmers separated by a newline (see also accompanying test files on github).
        If a list with only one file is given, returns the set of all kmers contained in that file
        If input files are a different format, this function can be adapted/overloaded
        """
        kmer_set = set()
        for file in file_list:
            with open(file, "r") as f:
                for line in f:
                    kmer_set.add(line.strip())
        return kmer_set

    def extract_quasi_primes(self, background_files: list, qp_files: list):
        logging.info(f"he following files are considered the target group and thus must contain all of the found quasi-primes: {qp_files}")
        logging.info(f"The following files are considered the background and thus cannot contain any of the found quasi-primes: {background_files}")

        # background_kmer_set is the union of all unique kmer peptides from all species except the species (or taxonomy) for which we want to extract
        # quasi-primes
        background_kmer_set = self.read_kmers(background_files)

        # observer_kmer_set is the union of all unique kmer peptides found in species (or taxonomy) for which we want to extract quasi-primes
        observed_kmer_set = self.read_kmers(qp_files)

        # kmers present in observed_kmer_set but not background_kmer_set are quasi primes for the observed species (or taxonomy)
        return observed_kmer_set.difference(background_kmer_set)

    def save_as_txt(self, kmer_set):
        with open(self.output_folder + self.timestamp + ".txt", "w") as f:
            for k in kmer_set:
                f.write(k + '\n')

    def calculate_intersections(self, kmer_files: list):
        """
        Calculate and store the intersection of all possible combinations of files given in the kmer_files list. Used for generating Venn diagrams for the
        accompanying paper. It is recommended to first combine kmers to a few big lists of respective categories (e.g. domains of life and viruses), as
        increasing the number of input files increases run time exponentially.
        """
        logging.info("Calculating the intersection of all possible combinations of the following kmer_files" + str(kmer_files))
        from itertools import combinations
        kmer_files_combinations = []
        for i in range(2, len(kmer_files) + 1):
            kmer_files_combinations += list(combinations(kmer_files, i))

        kmer_dict = {}
        for kmer_file in kmer_files:
            kmer_dict[kmer_file] = self.read_kmers([kmer_file])  # read_kmers expects a list

        with open(self.output_folder + self.timestamp + "_intersections.txt", "w") as f:
            for combination in kmer_files_combinations:
                # Calculate the intersection of all sets in # the combination
                intersection = kmer_dict[combination[0]].intersection(*[kmer_dict[x] for x in combination[1:]])
                f.write("Combination: " + str(combination) + "\n" + "Length of intersection: " + str(len(intersection)) + "\n" + str(intersection) + "\n")


if __name__ == "__main__":
    import glob
    BACKGROUND_FILES = glob.glob("test_files/*.5mers")
    QP_FILES = ["test_files/qp_test_file.5mersqp"]

    QPE = Quasi_Prime_Extractor("output/quasi_primes/")
    qp_set = QPE.extract_quasi_primes(BACKGROUND_FILES, QP_FILES)
    QPE.save_as_txt(qp_set)

    # Calculate and store intersections of all possible combinations of all kmer_files below
    QPE.calculate_intersections(["test_files/0_syntheticbacteria.5mers", "test_files/1_syntheticarchaea.5mers", "test_files/2_syntheticbacteria.5mers",
                                 "test_files/3_syntheticarchaea.5mers"])