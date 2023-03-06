import mmh3
import numpy as np

in_file = snakemake.input[0]
out_file = snakemake.output[0]


def hash_database(in_file, out_file, seed=18):
    """
    Perform mmh3 hashing on input database.
    MurmurHash3 is a non-cryptographic hashing algorithm.
    For more information on mmh3 algorithm visit https://pypi.org/project/mmh3/.

    :param in_file: str, input file
    :param out_file: str, output file
    :param seed: int, hashing seed: use identical integer for identical hashing results
    :return: database: np.array, hashed database
    """
    # initialize database
    database = np.empty([sum(1 for _ in open(in_file))])
    # hash protein accessions
    with open(in_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            database[line_num-1] = mmh3.hash64(line, signed=False, seed=seed)[0]
    # save in numpy format
    np.save(out_file, database)


def main():
    hash_database(in_file, out_file)


if __name__ == "__main__":
    main()