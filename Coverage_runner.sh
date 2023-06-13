MEMORY=100G
THREADS=10
wp=/media/HDD1/uri/projects/

bbmap.sh in1="$wp"Isra/Tuli/IS7_raw_assembly/IS7_combed_R1.fastq.gz \
        in2="$wp"Isra/Tuli/IS7_raw_assembly/IS7_combed_R2.fastq.gz \
        ref="$wp"Isra/Tuli/IS7_raw_assembly/48N_virus_scaffolds.fasta \
        covstats=48N_virus_es_scaffolds_covstats.txt  -Xmx"$MEMORY" threads=$THREADS  overwrite=t

bbmap.sh in1="$wp"Atlit_Isolates/48N/48N_trim_paired_R1.fastq \
        in2="$wp"Atlit_Isolates/48N/48N_trim_paired_R2.fastq \
        ref="$wp"Atlit_Isolates/48N/48N_hybrid_resolved.fasta \
        covstats=48N_old_reads_mapped_to_48N_hybrid_resolved.tsv  -Xmx"$MEMORY" threads=$THREADS  overwrite=t
