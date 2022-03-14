import matplotlib.pyplot as plt
import sys
import argparse

def parse_files(filenames):
	file_to_data = {}
	callers = set([])
	coverages = set([])
	for filename in filenames:
		precision = None
		recall = None
		fscore = None
		for line in open(filename, 'r'):
			if 'precision' in line:
				precision = float(line.split()[-1][:-1])
			if 'recall' in line:
				recall = float(line.split()[-1][:-1])
			if 'f1' in line:
				fscore = float(line.split()[-1][:-1])
		assert precision is not None
		assert recall is not None
		assert fscore is not None
		fields = filename.split('_')
		caller = fields[3]
		coverage = fields[2]
		file_to_data[(caller,coverage)] = [precision, recall, fscore]
		callers.add(caller)
		coverages.add(coverage)
	coverages = sorted([int(cov) for cov in list(coverages)])
	return file_to_data, sorted(list(callers)), [str(c) for c in coverages]


def plot_counts(stats, callers, coverages, outname, metric):
	plt.figure(figsize=(25,7))
	barWidth = 0.9

	m = None
	if metric == "precision":
		m = 0
	elif metric == "recall":
		m = 1
	else:
		m = 2

	n_callers = len(callers)
	n_coverages = len(coverages)
	for i,cov in enumerate(coverages):
		x_points = [1 + i + (n_coverages+1)*x  for x in range(n_callers)]
		y_points = [stats[(caller,cov)][m] for caller in callers ]
		print(x_points, y_points)
		plt.bar(x_points, y_points, width=0.9, label=cov)
	x_ticks = [0.5 + n_coverages/2 + i*(n_coverages+1) for i in range(n_callers)]
	plt.xticks(x_ticks, callers)
	txt_pos = [i + j*(n_coverages+1)  for j in range(n_callers) for i in range(1, 1+n_coverages)]
	vals = []
	for caller in callers:
		for cov in coverages:
			vals.append(stats[(caller,cov)][m])
	for pos, val in zip(txt_pos, vals):
		plt.text(x = pos -0.25, y=val, s=str(round(val, 2)), size=6)
	plt.ylabel(metric)
	plt.legend(loc=(1.04,0))
	plt.tight_layout()
	plt.savefig(outname)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='prec-rec_plots.py', description=__doc__)
	parser.add_argument('-c', '--callsets', nargs='+', required=True, help="truvari outputs for all combinations of callers/samples")
	parser.add_argument('-o', '--outname', required=True, help="name of the output file")
	parser.add_argument('-m', '--metric', required=True, choices=['precision', 'recall', 'fscore'])

	args = parser.parse_args()
	stats, callers, coverages = parse_files(args.callsets)
	print(stats)
	print(coverages)
	print(callers)

	plot_counts(stats, callers, coverages, args.outname, args.metric)
