import sys, argparse
import pandas as pd

parser = argparse.ArgumentParser(prog='filter-callset.py', description=__doc__)
parser.add_argument('-t', '--table', required=True, help='Table to which annotations shall be added.')
parser.add_argument('-o', '--outfile', required=True, help='The output table.')
parser.add_argument('-n', '--names', nargs='+', default=[], required=True, help='Names of datasets. Assumes cutesv>sniffles>svim.')
args = parser.parse_args()


df = pd.read_csv(args.table, sep='\t')

# filtered set: remove variants only seen in svim/sniffles and keep rest
for name in args.names:
	colname = name + '_filtered'
	df[colname] = ((1*df['in_cutesv_' + name]) + (1*df['in_sniffles_' + name]) + (1*df['in_svim_' + name]) ) >= 2
	nr_filtered = df[colname].sum()
	print('Filtered set for ' + name + ' contains ' + str(nr_filtered) + ' variants.')
	colname_union = name + '_union'
	df[colname_union] = df['in_cutesv_' + name] | df['in_sniffles_' + name] | df['in_svim_' + name]
	nr_union = df[colname_union].sum()
	print('Union set for ' + name + ' contains ' + str(nr_union) + ' variants.')
df.to_csv(args.outfile, sep='\t', index=False, na_rep='nan')
