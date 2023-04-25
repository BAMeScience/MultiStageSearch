import pandas as pd
import matplotlib.pyplot as plt


def CreateStrainBarPlots(in_file, out_file):
    df = pd.read_csv(in_file, delimiter="\t")

    df["strain"] = df["strain"].astype(str)
    counts = df["counts"]
    strain = df["strain"]

    fig, ax = plt.subplots()
    ax.bar(strain, counts)

    # Add labels and title
    ax.set_xlabel("Strain")
    ax.set_ylabel("Counts")
    ax.set_title("Bar Plot of Counts vs Strain")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateTaxIdScoresBarPlot(in_file, out_file):
    df = pd.read_csv(in_file, sep="\t")

    taxids = df["taxid"].astype(str)[:20]
    weights = df["weight"][:20]

    plt.bar(taxids, weights)
    plt.xlabel("Tax ID")
    plt.ylabel("Weight")
    plt.title("Tax ID vs Weight")

    plt.xticks(rotation=90)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateProportionsPieChart(first_search, final_search, out_file):
    df1 = pd.read_csv(first_search, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
    df1.columns = ["Peptide"]
    df1["Peptide"] = df1.Peptide.apply(lambda x: x.split(","))
    report_df_1 = df1.explode("Peptide", ignore_index=True)

    df2 = pd.read_csv(final_search, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
    df2.columns = ["ORF"]
    df2["ORF"] = df2.ORF.apply(lambda x: x.split(","))
    report_df_2 = df2.explode("ORF", ignore_index=True)

    num_entries = [len(report_df_1), len(report_df_2)]
    labels = ['First Search', 'Final Search']

    plt.pie(num_entries, labels=labels, autopct='%1.1f%%')
    plt.title('Number of PSMs')
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def CreateConfidenceHistogram(in_file, out_file):
    df = pd.read_csv(in_file, sep="\t", on_bad_lines="skip")
    df.rename(columns={"Protein(s)": "Proteins"}, inplace=True)
    df["Proteins"] = df.Proteins.apply(lambda x: x.split(","))
    exploded_df = df.explode("Proteins", ignore_index=True)

    plt.hist(exploded_df['Confidence [%]'], bins=20)
    plt.xlabel('Confidence')
    plt.ylabel('Count')
    plt.title('Confidence Histogram')
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.clf()


def main():
    strain_mappings = snakemake.input[0]
    taxID_scores = snakemake.input[1]
    frist_search = snakemake.input[2]
    final_search = snakemake.input[3]

    strain_bar_plot = snakemake.output[0]
    taxIdScores_bar_plot = snakemake.output[1]
    proportions_pie_chart = snakemake.output[2]
    first_search_confidence_histogram = snakemake.output[3]
    final_search_confidence_histogram = snakemake.output[4]

    CreateStrainBarPlots(strain_mappings, strain_bar_plot)
    CreateTaxIdScoresBarPlot(taxID_scores, taxIdScores_bar_plot)
    CreateProportionsPieChart(frist_search, final_search, proportions_pie_chart)
    CreateConfidenceHistogram(frist_search, first_search_confidence_histogram)
    CreateConfidenceHistogram(final_search, final_search_confidence_histogram)


if __name__ == "__main__":
    main()