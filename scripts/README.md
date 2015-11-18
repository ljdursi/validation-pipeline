Directory for the scripts which make up the actual nuts and bolts of
the pipeline.  Most of the scripts in here only operate on one item
or matched pair, with the `run-*` scripts above launching them over
all available samples

* `bamToFastq.sh`, `align_one_fastq.sh`, `alignfastqs.sh` - ingest BAM, convert to FASTQ, and realign with pancancer standard tools, reference
* `bam-depth.sh`, `bam-stats.sh`, `runfastqc.sh` - QC stats on realigned BAMs
* `build-beds.sh`, `union-germline.sh` - take selected or germline calls and build BED files surounding them, or merge germline indel calls from multiple callers
* `one-indelrealign.sh` - realign remapped BAM file around possible germline indels
* `generate-IGV-scripts.sh` - generate the batch scripts for IGV to visualize the bams around each selected variant
* `one-readcount.sh`+`snv_readcounts.py`, `one-sga-annotate.sh`, `count-sv-support.py` - get counts of support for and against purported variants for a sample, for SNVs, Indels, and SVs, respectively
* `one-call.sh`+`snv_indel_call.py` - make calls for SNVs or Indels based on the count support above
* `add-germline-indel-dist-info.sh` - add germline nearest-indel information to a VCF
* `merge-validation-selection.sh` - merge the validation calls VCF with the master VCFS (with WGS context information) 
* `parsecombined_snv.py` - convert the above merged VCFs into a CSV suitable for analysis
