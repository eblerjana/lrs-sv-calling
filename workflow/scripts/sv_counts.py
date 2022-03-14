import matplotlib.pyplot as plt
import sys



###########################
# Plot Sv counts for all samples sorted by SV caller
#
# Call:
#  python code/sv_counts.py [-log or -nolog] [ [filename_sample1, filename_sample2, ...] for all sv_callers]  output_file
#
# Handles variable number of and samples
# (samples between sv_callers must be same and in same order)
# Only hardcoded colors for up to 5 samples. add more if you need more
###########################


def get_counts(fnames):
    counts = [0 for _ in fnames]
    for i, fname in enumerate(fnames):

        print("Counting file", fname)
        with open(fname,"r") as infile:

            # count svs of single file
            for line in infile:
                # header
                if line.startswith("#"):
                    continue

                # count
                infos = line.split()[7]
                if "LEN=" in infos:
                    svlen = infos.split("LEN=")[1].split(";")[0]
                    if abs(int(svlen)) >= 50:
                        counts[i] += 1
                elif "END=" in infos:
                    svlen = int(infos.split("END=")[1].split(";")[0]) - int(line.split()[1])
                    if svlen >= 50:
                        counts[i] += 1
    return counts


# sv caller order follows input order
# sample order follows input order
def plot_counts(all_sv_counts, num_svcaller, num_samples, sample_names, caller_names, log, out_fname):
    plt.figure(figsize=(25,7))
    barWidth = 0.9

    if num_samples <= 5:
        colors = ["darkblue", "mediumblue", "blue", "royalblue", "cornflowerblue"]

    # plot bars by sample
    for i in range(num_samples):
        x_points = [1 + i + (num_samples+1)*x  for x in range(num_svcaller)]
        y_heights = [y[i] for y in all_sv_counts]
        if num_samples <= 5:
            plt.bar(x_points, y_heights, color = colors[i], width=barWidth, label = sample_names[i])
        else:
            plt.bar(x_points, y_heights, width=barWidth, label = sample_names[i])

    # plot x ticks for every sv caller
    x_ticks = [0.5 + num_samples/2 + i*(num_samples+1)  for i in range(num_svcaller)]
    print(x_ticks)
    plt.xticks(x_ticks, caller_names)

    # add values on top of bars
    txt_pos = [i + j*(num_samples+1)  for j in range(num_svcaller) for i in range(1, 1+num_samples)]
    print(txt_pos)
    vals = [count for sv_counts in all_sv_counts for count in sv_counts]
    for pos, val in zip(txt_pos, vals):
        plt.text(x = pos -0.25 , y = val , s = str(val), size = 6)

    # y scale logarithmic
    plt.legend(loc=(1.04,0))
    if log:
        plt.yscale("log")
        plt.ylabel("SV count (log)")
    else:
        plt.ylabel("SV count")
    plt.tight_layout()


    # save
    print("Plotting file", out_fname)
    plt.savefig(out_fname)


def parse_argv(argv):
    out_fname = argv[-1]
    log = argv[1] == "-log"
    argv = argv[2:-1]

    num_svcaller = 3
    num_samples = len(argv)//3

    sample_names = [fname.split("/")[-1].split('_')[2] for fname in argv[: num_samples]]
    caller_names = [argv[i * num_samples].split("/")[-1].split('_')[-2] for i in range(num_svcaller)]

    print("Samples:", sample_names)
    print("Number of samples:", num_samples)
    print("Number of SV callers:", num_svcaller)

    return out_fname, log, argv, num_svcaller, num_samples, sample_names, caller_names


def main():
    out_fname, log, argv, num_svcaller, num_samples, sample_names, caller_names = parse_argv(sys.argv)

    # Count all svs of all files
    # all_sv_counts = [ [count_sample1, count_sample2, ....]  for each sv caller]
    all_sv_counts = []
    for i in range(num_svcaller):
        all_sv_counts.append(get_counts(argv[i * num_samples: (i+1) * num_samples]))

    print("All SV counts:")
    print(all_sv_counts)

    plot_counts(all_sv_counts, num_svcaller, num_samples, sample_names, caller_names, log, out_fname)



if __name__ == "__main__":
    main()
