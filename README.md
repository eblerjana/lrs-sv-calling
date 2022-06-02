# SV calling pipeline using long sequencing reads

The pipeline calls structural variants from aligned long sequencing reads (PacBio HiFi, ONT, ONT ultralong) using three SV callers: cuteSV, sniffles and SVIM. It additionally filters the calls, creates a merged callset and generates some plots comparing the sets. Variant calls for HG002 are additionally evaluated based on the GIAB medically relevant SVs.

## How to run

For each BAM file to be analyzed, the path to the corresponding BAM file must be specified in config/config.yaml, together with a sample name, the coverage and the sequencing technology ("hifi", "ont" or "ont-ul").
Furthermore, the path to the reference sequence used must be added to the config file.
