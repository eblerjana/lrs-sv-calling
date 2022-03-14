import argparse
import sys

"""
usage: filter_vcf_qual_size.py [-h] [-q Q] input_file output_file

positional arguments:
  input_file   The sv file (.vcf)
  output_file  The filtered sv file (.vcf)

optional arguments:
  -h, --help   show this help message and exit
  -q Q         Quality cutoff for filtering as integer
"""


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="The sv file (.vcf)")
    parser.add_argument("output_file", help="The filtered sv file (.vcf)")
    parser.add_argument("-q", help="Quality cutoff for filtering as integer", type=int, default=10)
    return parser.parse_args()


def filter_file(infname, q, outfname):
    with open(infname,"r") as infile, open(outfname,"w") as outfile:
        print("Filtering file", infname)
        line_count = 0
        keep_count = 0
        
        drop_invalid_start = 0
        drop_no_length = 0
        drop_low_qual = 0
        drop_types = 0
        drop_short = 0
        drop_long = 0

        # filter line by line
        for line in infile:

            # keep header
            if line.startswith("#"):
                outfile.write(line)
                continue

            # count non header lines (SVs)
            line_count +=1
            columns = line.split()
            
            # skip if quality is given and quality < q
            if columns[5] == ".":
                pass
            elif float(columns[5]) < q:
                drop_low_qual += 1
                continue

            info_field = {i.split('=')[0]:i.split('=')[1] for i in columns[7].split(';') if '=' in i if '=' in i}
            
            # only keep specific variant types
            assert "SVTYPE" in info_field
            if not info_field["SVTYPE"] in ["INS", "DEL", "INV", "DUP", "INVDUP", "CNV", "DUP:TANDEM"]:
            	sys.stderr.write("Found " + info_field["SVTYPE"] + "\n")
            	drop_types += 1
            	continue

            # check if start position is valid (VCF is 1-based)
            if columns[1] == "0":
                print('INVALID:', line)
                drop_invalid_start += 1
                continue

            
            # skip if variant does not have a length or end
            if not "SVLEN" in info_field and not "END" in info_field:
                drop_no_length += 1
                continue

            # skip if length is < 50 or > 500000
            if "SVLEN" in info_field:
                svlen = abs(int(info_field['SVLEN']))
                if svlen < 50:
                    drop_short += 1
                    continue
                if svlen > 500000:
                    drop_long += 1
                    continue
            else:
                svlen = int(info_field['END']) - int(columns[1])
                assert svlen > 0
                if svlen < 50:
                    drop_short += 1
                    continue
                if svlen > 500000:
                    drop_long += 1
                    continue

            # if line not dropped by previous checks, keep it
            outfile.write(line)
            keep_count +=1

        assert drop_invalid_start + drop_low_qual + drop_types + drop_short + drop_long + drop_no_length + keep_count == line_count 
        print("Total SV count:", line_count)
        print("Dropped due to invalid start coordinate:", drop_invalid_start)
        print("Dropped due to low quality (<" + str(q) + "):", drop_low_qual)
        print("Dropped due to invalid type:", drop_types)
        print("Dropped due to missing length/coordinates:", drop_no_length)
        print("Dropped due to short length (<50):", drop_short)
        print("Dropped due to long length (>500000):", drop_long)
        print("Kept SVs:", keep_count)



def main():
    args = parse_arguments()
    filter_file(args.input_file, args.q, args.output_file)


if __name__ == "__main__":
    main()
