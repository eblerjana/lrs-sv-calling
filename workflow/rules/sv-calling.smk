import math
configfile: "config/config.yaml"


### the input alignments are expected in: data/hifi/ and data/ont/ ###

info_to_filename = {}
samples = set([])
tech = set([])
coverages = set([])

for filename in config['data']:
	info = [i for i in config['data'][filename]]
	sample = info[0]
	coverage = info[1]
	t = info[2]
	samples.add(sample)
	assert t in ['hifi', 'ont', 'ont-ul']
	tech.add(t)
	coverages.add(int(coverage))
	info_to_filename[(sample, str(coverage), t)] = filename

samples = list(samples)
tech = list(tech)
coverages = [str(c) for c in sorted(list(coverages))]

print(samples, tech, coverages)
print(info_to_filename)


# cuteSV needs different parameters for hifi/ont data
cutesv_params = {
	"hifi": "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5",
	"ont": "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3",
	"ont-ul": "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3"
}


# determine cuteSV min_support cutoffs and quality cutoffs for filtering svim calls
cutesv_min_support = {}
cutoff_qual = {}
#for cov in coverages:
#	cutesv_min_support[cov] = int(min(10, round(int(cov)*0.6, 0)))  #int(min(10, math.ceil(int(cov)*0.25)))
#	cutoff_qual[cov] = math.ceil(int(cov)/4.0)

for cov in coverages:
	if int(cov) <= 5:
		cutesv_min_support[cov] = 2
		cutoff_qual[cov] = 2
	elif int(cov) <= 10:
		cutesv_min_support[cov] = 3
		cutoff_qual[cov] = 3
	elif int(cov) <= 20:
		cutesv_min_support[cov] = 4
		cutoff_qual[cov] = 4
	elif int(cov) <= 30:
		cutesv_min_support[cov] = 5
		cutoff_qual[cov] = 5
	else:
		cutesv_min_support[cov] = 10
		cutoff_qual[cov] = 10



############################
# download GIAB CMRG data
############################


rule download_giab_med_svs:
	output:
		vcf="results/data/giab-med/giab-med.vcf.gz",
		tbi="results/data/giab-med/giab-med.gz.tbi",
		bed="results/data/giab-med/giab-med.bed"
	shell:
		"""
		wget -O {output.vcf} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
		wget -O {output.tbi} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi
		wget -O {output.bed} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed
		"""



#########################
#     Sniffles rule     #
#########################

rule sniffles:
	input:
		bam= lambda wildcards: info_to_filename[(wildcards.sample, wildcards.coverage, wildcards.tech)],
		reference = config['reference']
	output:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_sniffles.vcf"
	conda:
		"../envs/sniffles_env.yaml"
	wildcard_constraints:
		sample = "|".join(samples),
		tech = "|".join(tech)
	log:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_sniffles.log"
	benchmark:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_sniffles_benchmark.txt"
	threads: 24
	resources:
		mem_total_mb=50000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		"sniffles -t {threads} -i {input.bam} -v {output} --reference {input.reference} --minsvlen 50 &> {log}"


#########################
#       svim rule       #
#########################


# call svs with svim
rule svim:
	input:
		bam= lambda wildcards: info_to_filename[(wildcards.sample, wildcards.coverage, wildcards.tech)],
		reference=config["reference"]
	output:
		raw="results/{tech}/{sample}/{tech}_{sample}_{coverage}/variants.vcf",
		vcf="results/{tech}/{sample}/{tech}_{sample}_{coverage}_svim.vcf"
	conda:
		"../envs/svim_env.yaml"
	wildcard_constraints:
		sample = "|".join(samples),
		tech = "|".join(tech)
	log:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_svim.log"
	benchmark:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_svim_benchmark.txt"
	threads: 1
	resources:
		mem_total_mb=20000,
		runtime_hrs=8,
		runtime_min=1
	params:
		outdir = "results/{tech}/{sample}/{tech}_{sample}_{coverage}/"
	shell:
		"""
		svim alignment {params.outdir} {input.bam} {input.reference} --min_sv_size 50 &>> {log}
		cp {output.raw} {output.vcf}
		"""



#########################
#     CuteSV rule       #
#########################

rule cutesv:
	"""
	NB: cuteSV fails if temp working directory does not exist at execution time
	"""
	input:
		bam= lambda wildcards: info_to_filename[(wildcards.sample, wildcards.coverage, wildcards.tech)],
		reference=config["reference"]
	output:
		vcf="results/{tech}/{sample}/{tech}_{sample}_{coverage}_cutesv.vcf",
		tmp=temp(directory("results/{tech}/{sample}/{tech}_{sample}_{coverage}_cutesv.tmp.wd/"))
	log:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_cutesv.log"
	benchmark:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_cutesv_benchmark.txt"
	conda:
		"../envs/cutesv_env.yaml"
	wildcard_constraints:
		sample = "|".join(samples),
		tech = "|".join(tech)
	threads: 24
	resources:
		mem_total_mb=20000,
		runtime_hrs=8,
		runtime_min=1
	params:
		out_dir = lambda wildcards, output: output.vcf.replace('.vcf', '.tmp.wd'),
		call_params = lambda wildcards: cutesv_params[wildcards.tech],
		min_supp = lambda wildcards: cutesv_min_support[wildcards.coverage]
	shell:
		"""
		mkdir -p {params.out_dir}
		cuteSV -t {threads} -S {wildcards.sample} {params.call_params} {input.bam} {input.reference} {output.vcf} --genotype -l 50 -s {params.min_supp} {params.out_dir} &> {log}
		"""


##################################
#  filter out low quality calls  #
##################################


rule filter_calls:
	"""
	filter out low quality calls. This is necessary for svim, because
	it reports all calls (including low quality ones). coverage / 4.0 is
	the author's suggestion.
	"""
	input:
		vcf = "results/{tech}/{sample}/{tech}_{sample}_{coverage}_{caller}.vcf"
	output:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_{caller}_filtered.vcf"
	log:
		"results/{tech}/{sample}/{tech}_{sample}_{coverage}_{caller}_filtered.log"
	wildcard_constraints:
		samples = "|".join(samples),
		tech = "|".join(tech),
		caller = "cutesv|sniffles|svim"
	params:
		cutoff = lambda wildcards: cutoff_qual[wildcards.coverage] if wildcards.caller == "svim" else 0
	threads: 1
	shell:
		"python workflow/scripts/filter_vcf_qual_size.py -q {params.cutoff} {input.vcf} {output} &>> {log}"



####################################
# compute precision/recall
####################################

rule compress_vcf:
	input:
		"{filename}.vcf"
	output:
		"{filename}.vcf.gz"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""


rule precision_recall:
	input:
		vcf="results/{tech}/{sample}/{tech}_{sample}_{coverage}_{caller}_filtered.vcf.gz",
		truth="results/data/giab-med/giab-med.vcf.gz",
		bed="results/data/giab-med/giab-med.bed",
		reference=config["reference"]
	output:
		"results/{tech}/{sample}/evaluation/{tech}_{sample}_{coverage}_{caller}_truvari/summary.txt"
	params:
		tmp="results/{tech}/{sample}/evaluation/{tech}_{sample}_{coverage}_{caller}_truvari_temp",
		outname="results/{tech}/{sample}/evaluation/{tech}_{sample}_{coverage}_{caller}_truvari",
		p= lambda wildcards: "--pctsim=0" if wildcards.caller == "svim" else "" # svim does not generate sequence resolved calls
	resources:
		mem_total_mb = 30000,
	log:
		"results/{tech}/{sample}/evaluation/{tech}_{sample}_{coverage}_{caller}_truvari.log"
	shell:
		"""
		truvari bench -b {input.truth} -c {input.vcf} -f {input.reference} -o {params.tmp} --multimatch -r 2000 --no-ref a -C 2000 --includebed {input.bed} {params.p} --passonly &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""



####################################
#  compare callsets
####################################


def merge_callsets_files(wildcards):
	result = []
	for coverage in coverages[::-1]:
		for caller in ["cutesv", "sniffles", "svim"]:
			result.append("results/{tech}/{sample}/{tech}_{sample}_{coverage}_{caller}_filtered.vcf".format(tech=wildcards.tech, sample=wildcards.sample, coverage=coverage, caller=caller))
	return result


def merge_callsets_names(wildcards):
	result = []
	for coverage in coverages[::-1]:
		for caller in ["cutesv", "sniffles", "svim"]:
			result.append("_".join([caller, wildcards.sample, coverage]))
	return result


# compute intersection between callers/samples
rule merge_callsets:
	input:
		merge_callsets_files
	output:
		tsv = temp("results/{tech}/{sample}/intersection/{tech}_{sample}_intersection_raw.tsv"),
		vcf = temp("results/{tech}/{sample}/intersection/{tech}_{sample}_intersection_raw.vcf"),
		pdf= "results/{tech}/{sample}/intersection/{tech}_{sample}_intersection.pdf"
	log:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_intersection.log"
	benchmark:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_intersection_bechmark.txt"
	params:
		names =  merge_callsets_names
	wildcard_constraints:
		tech = "hifi|ont|ont-ul"
	resources:
		mem_total_mb=10000,
		runtime_hrs=20,
		runtime_min=1
	conda:
		"../envs/plot_env.yaml"
	shell:
		"python3 workflow/scripts/intersect_callsets.py intersect -t {output.tsv} -v {output.vcf} -p {output.pdf} -c {input} -n {params.names} &> {log}"


# sort callsets
rule sort_merged_callset:
    input:
        tsv="results/{tech}/{sample}/intersection/{tech}_{sample}_intersection_raw.tsv",
        vcf="results/{tech}/{sample}/intersection/{tech}_{sample}_intersection_raw.vcf"
    output:
        tsv="results/{tech}/{sample}/intersection/{tech}_{sample}_intersection.tsv",
        vcf="results/{tech}/{sample}/intersection/{tech}_{sample}_intersection.vcf"
    shell:
    	"""
    	cat {input.tsv} | awk '$1 ~ /^ID/ {{print $0;next}} {{print $0 | "sort -k2,2 -k3,3n"}}' > {output.tsv}
    	cat {input.vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' > {output.vcf}
    	"""



############################################
#  integrate callsets to define final set
############################################


# define filtered callsets for each coverage that contains all variants calle by at least two callers
rule filtered_callset:
	input:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_intersection.tsv"
	output:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_integrated.tsv",
	conda:
		"../envs/plot_env.yaml"
	log:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_integrated.log",
	params:
		names = lambda wildcards: ["_".join([wildcards.sample, cov]) for cov in coverages]
	shell:
		"python3 workflow/scripts/filter_callsets.py -t {input} -n {params.names} -o {output}"



############################################
#  create plots
############################################


# plots number of SV counts per sample and caller
rule plot_sv_counts:
	input:
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_cutesv_filtered.vcf", coverage=coverages),
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_sniffles_filtered.vcf", coverage=coverages),
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_svim_filtered.vcf", coverage=coverages)
	output:
		"results/{tech}/{sample}/plots/{tech}_{sample}_sv-counts.pdf"
	conda:
		"../envs/plot_env.yaml"
	log:
		"results/{tech}/{sample}/plots/{tech}_{sample}.log"
	threads: 1
	shell:
		"python workflow/scripts/sv_counts.py -nolog {input} {output} &>> {log}"


# plots length distribution of SV callsets
rule plot_len_dist:
	input:
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_cutesv_filtered.vcf", coverage=coverages),
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_sniffles_filtered.vcf", coverage=coverages),
		expand("results/{{tech}}/{{sample}}/{{tech}}_{{sample}}_{coverage}_svim_filtered.vcf", coverage=coverages)
	output:
		"results/{tech}/{sample}/plots/{tech}_{sample}_sv-lengths.pdf"
	conda:
		"../envs/plot_env.yaml"
	log:
		"results/{tech}/{sample}/plots/{tech}_{sample}_sv-lengths.log"
	threads: 1
	shell:
		"python workflow/scripts/sv_len_distr.py {input} {output} &>> {log}"


# plot evaluation results
rule plot_prec_rec:
	input:
		expand("results/{{tech}}/{{sample}}/evaluation/{{tech}}_{{sample}}_{coverage}_{caller}_truvari/summary.txt", coverage=coverages, caller=["cutesv", "sniffles", "svim"])
	output:
		"results/{tech}/{sample}/plots/{tech}_{sample}_{metric}.pdf"
	wildcard_constraints:
		metric="precision|recall|fscore"
	conda:
		"../envs/plot_env.yaml"
	shell:
		"python3 workflow/scripts/prec-rec_plots.py -c {input} -o {output} --metric {wildcards.metric}"


# for a caller, upset plot across all coverages
rule plot_upset_across_files:
	input:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_integrated.tsv",
	output:
		"results/{tech}/{sample}/plots/{tech}_{sample}_{caller}_upset.pdf"
	wildcard_constraints:
		caller = "cutesv|sniffles|svim" 
	conda:
		"../envs/plot_env.yaml"
	params:
		columns = lambda wildcards: ["_".join(["in", wildcards.caller, wildcards.sample, cov]) for cov in coverages]
	shell:
		"python3 workflow/scripts/plot_upset.py -t {input} -o {output} -n {params.columns}"


# upset plot for filtered set
rule plot_upset_filtered_callsets:
	input:
		"results/{tech}/{sample}/intersection/{tech}_{sample}_integrated.tsv",
	output:
		"results/{tech}/{sample}/plots/{tech}_{sample}_upset-integrated.pdf"
	conda:
		"../envs/plot_env.yaml"
	params:
		columns = lambda wildcards: ["_".join([wildcards.sample, cov, "filtered"]) for cov in coverages]
	shell:
		"python3 workflow/scripts/plot_upset.py -t {input} -o {output} -n {params.columns}"
