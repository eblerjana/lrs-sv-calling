import matplotlib.pyplot as plt
import numpy as np
import sys


###########################
# Plot SV length distribution for all given files
#
# Call:
#  python code/sv_len_distr.py [ [filename_sample1, filename_sample2, ...] for all sv_callers]  output_file
#
# (samples between sv_callers must be same and in same order)
# Hardcoded colors. add more if you need more
# also might need to adjust the lables in the plot function
###########################



def get_lengths(fname):
    print("Counting file", fname)
    lengths = []
    with open(fname,"r") as infile:
        for line in infile:

            #header
            if line.startswith("#"):
                continue

            # find length
            infos = line.split()[7]
            if "LEN=" in infos:
                svlen = infos.split("LEN=")[1].split(";")[0]
                svlen = abs(int(svlen))  # take the absolute value
                if svlen < 50:
                    continue
                lengths.append(svlen)
            elif "END=" in infos:
                svlen = int(infos.split("END=")[1].split(";")[0]) - int(line.split()[1])
                if svlen < 50:
                    continue
                lengths.append(svlen)

    return lengths



def plot_dists(all_lengths, num_svcaller, num_samples, names, out_fname):
    plt.figure(figsize=(25,7))

    # lables only work for the right sample names
    labels = names

    blue = ["darkblue", "mediumblue", "blue", "royalblue", "cornflowerblue", "slateblue", "darkslateblue", "deepskyblue", "dodgerblue", "steelblue"]
    green = ["darkgreen","forestgreen","green", "limegreen", "lime", "mediumseagreen", "mediumspringgreen", "aquamarine", "greenyellow", "lightgreen"]
    red = ["darkred", "brown", "indianred", "salmon", "lightsalmon", "lightcoral", "firebrick", "orangered", "red", "chocolate"]
    palette = [blue, green, red]

    colors = []
    assert num_svcaller < 4
    for i in range(num_svcaller):
        colors += palette[i][0:num_samples]

    # plot distributions with logarithmic bins on x axis
    for lengths, label, color in zip(all_lengths, labels, colors):
        logbins = np.geomspace(min(lengths), max(lengths), 70)
        plt.hist(lengths, histtype="step",bins =  logbins,  label=label, color = color)
        #plt.hist(lengths, histtype="step",bins =  70,  label=label, color = color) # for non log x axis

    # log scales
    plt.xscale("log")
    plt.yscale("log")

    # save fig
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.legend(loc=(1.04,0))
    plt.tight_layout()
    plt.savefig(out_fname)



def parse_argv(argv):
    out_fname = argv[-1]
    argv = argv[1:-1]

    num_svcaller = 3
    num_samples = len(argv)//3

#    sample_names = [fname.split("/")[-1][:-4] for fname in argv[: num_samples]]
#    caller_names = [argv[i * num_samples].split("/")[1] for i in range(num_svcaller)]

    names = [fname.split('/')[-1][:-4] for fname in argv]
    print("Number of samples:", num_samples)
    print("Number of SV callers:", num_svcaller)

    return out_fname, argv, num_svcaller, num_samples, names



def main():
    out_fname, argv, num_svcaller, num_samples, names = parse_argv(sys.argv)

    # collect all lengths
    # all_lengths = [ [lengths of file1], [lengths of file2], ....  for each input file]
    all_lengths = []
    for fname in argv:
        all_lengths.append(get_lengths(fname))

    plot_dists(all_lengths, num_svcaller, num_samples, names, out_fname)




if __name__ == "__main__":
    main()
