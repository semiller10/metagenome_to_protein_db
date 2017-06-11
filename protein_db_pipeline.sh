#!/bin/bash

# List the metagenome urls
urls=(
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1019366/ERR1019366_1.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1019366/ERR1019366_2.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/ERA486/ERA486537/fastq/UOA_33_C6-Cheng_TTAGGC_L005_R1.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/ERA486/ERA486537/fastq/UOA_33_C6-Cheng_TTAGGC_L005_R2.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR102/007/ERR1022687/ERR1022687_1.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR102/007/ERR1022687/ERR1022687_2.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/004/ERR1034454/ERR1034454_1.fastq.gz \
    ftp.sra.ebi.ac.uk/vol1/fastq/ERR103/004/ERR1034454/ERR1034454_2.fastq.gz
    )
mgf_files=(
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_2_2_1_1_sem.mgf \
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_2_4_1_1_sem.mgf \
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_26_2_1_1_sem.mgf \
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_26_4_1_1_sem.mgf \
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_27_2_1_1_sem.mgf \
    /home/samuelmiller/toolik_proteomes/042017_toolik_core_27_4_1_1_sem.mgf
    )
cores=24
python_module=python/3.5-2017q2
working_dir="$(pwd)"
mkdir -p trimmed translated
SolexaQA_min_nt_prob=0.01
SolexaQA_bin="$(find ~ -type f -name SolexaQA++)"
megahit_bin="$(find ~ -type f -name megahit)"
megahit_out_dir="$working_dir"/megahit_out
megahit_k_min=27
megahit_k_max=87
megahit_k_step=10
assembly="$megahit_out_dir"/intermediate_contigs/k"$megahit_k_max".contigs.fna
megahit_toolkit_bin="$(find ~ -type f -name megahit_toolkit)"
Graph2Pro_dir="$(find ~ -type d -name Graph2Pro)"
remove_500bp_script="$Graph2Pro_dir"/utils/remove_500bp.py
FGSP_dir="$(find ~ -type d -name FragGeneScanPlus)"
is_whole_gene=0
FGSP_model=illumina_10
FGSP_alloc=3000
FGSP_output_metadata=1
DBGraph2Pro_bin="$Graph2Pro_dir"/DBGraph/DBGraph2Pro
DBGraph2Pro_depth=5
createFixedReverseKR_script="$Graph2Pro_dir"/utils/createFixedReverseKR.py
MSGFP_jar="$Graph2Pro_dir"/MSGF+/MSGFPlus.jar
java_alloc=-Xmx32g
instrument_type=0
parent_mass_tol=10ppm
isotope_error_range=-1,2
mods_file="$Graph2Pro_dir"/MSGF+/standard_mods.txt
is_target_decoy_db=0
min_precursor_charge=1
max_precursor_charge=4
MSGFP_mzid_to_tsv=edu.ucsd.msjava.ui.MzIDToTsv
report_decoy_psms=1
parseFDR_script="$Graph2Pro_dir"/utils/parseFDR.py
fdr_level=0.05
combineFragandDBGraph_script="$Graph2Pro_dir"/utils/combineFragandDBGraph.py
DBGraphPep2Pro_bin="$Graph2Pro_dir"/DBGraph/DBGraphPep2Pro

# Loop through each metagenome
for (( i=0; i<"${#urls[@]}"; i+=2 )) ; do
    
    paired_fasta_files=()
    # Loop through each of 2 fastq files of paired-end reads
    for j in "$(seq "$i" "$(($i+1))")" ; do
        # Download fastq files
        url="${urls[j]}"
        url_array=( ${url//\// } )
        fastq_filename="${url_array[-1]}"
        wget -O "$fastq_filename" "$url"

        # Quality control
        "$SolexaQA_bin" dynamictrim "$fastq_filename" -d trimmed -p "$SolexaQA_min_nt_prob"
        fastq_filename_array=( ${fastq_filename//\./ } )
        fastq_basename="${fastq_filename_array[0]}"
        sed -n '1~4s/^@/>/p; 2~4p' trimmed/"$fastq_filename".trimmed > trimmed/"$fastq_basename".trimmed.fna

        paired_fasta_files+=( trimmed/"$fastq_basename".trimmed.fna )

        # Delete downloaded and trimmed fastq files
        rm "$fastq_filename"
        rm trimmed/"$fastq_filename".trimmed
    done

    # Make assembly
    "$megahit_bin" -1 "${paired_fasta_files[0]}" -2 "${paired_fasta_files[1]}" \
    --k_min "$megahit_k_min" --k_max "$megahit_k_max" --k_step "$megahit_k_step" \
    -t "$cores" -o "$megahit_out_dir"

    # Delete all assembly output except the last intermediate contigs
    mv "$assembly" "$working_dir"
    assembly="$working_dir"/k"$megahit_k_max".contigs.fna
    rm -r "$megahit_out_dir"

    # Make assembly graph
    metagenome_name="${fastq_basename%_2}"
    assembly_graph="$working_dir"/"$metagenome_name"."$megahit_k_max".fastg
    "$megahit_toolkit_bin" contig2fastg "$megahit_k_max" \
    "$assembly" > "$assembly_graph"

    # Graph2Pro
    # 1. Remove contigs shorter than 500 nt
    module load $python_module
    long_assembly="$metagenome_name"."$megahit_k_max".500.contigs.fna
    python "$remove_500bp_script" "$assembly" > "$long_assembly"
    rm "$assembly"
    
    # 2. Translate contigs with FragGeneScan+
    cd "$FGSP_dir"
    long_contig_peps="$working_dir"/"$metagenome_name"."$megahit_k_max".500.faa
    long_contig_pep_dna="$working_dir"/"$metagenome_name"."$megahit_k_max".500.ffn
    ./FGS+ -s "$long_assembly" -o "$working_dir"/"$metagenome_name"."$megahit_k_max".500 \
    -w "$is_whole_gene" -t "$FGSP_model" -p "$cores" -m "$FGSP_alloc" -d "$FGSP_output_metadata"

    # 3. Find all tryptic peptides in assembly graph with Graph2Pep
    cd "$working_dir"
    assembly_graph_peps="$working_dir"/"$metagenome_name".contig_pep.fasta
    "$DBGraph2Pro_bin" -d "$DBGraph2Pro_depth" \
    -e "$assembly_graph" -s "$assembly" -o "$assembly_graph_peps"

    # 4. Create decoy dbs from peptides
    long_contig_peps_decoy_db="$working_dir"/"$metagenome_name".long_contig_peps.fixedKR.fasta
    python "$createFixedReverseKR_script" "$long_contig_peps_decoy_db"

    assembly_graph_peps_decoy_db="$working_dir"/"$metagenome_name".assembly_graph_peps.fixedKR.fasta
    python "$createFixedReverseKR_script" "$assembly_graph_peps_decoy_db"

    # Loop through each proteome
    for mgf_file in ${mgf_files[@]} ; do

        # 5.a. Spectra search against peptide databases
        long_contig_peps_mzid_psms="$working_dir"/"$metagenome_name".long_contig_peps.fixedKR.mzid
        java "$java_alloc" -jar "$MSGFP_jar" -s "$mgf_file" -d "$long_contig_peps_decoy_db" -o "$long_contig_peps_mzid_psms" \
        -inst "$instrument_type" -t "$parent_mass_tol" -ti "$isotope_error_range" \
        -mod "$mods_file" -tda "$is_target_decoy_db" \
        -minCharge "$min_precursor_charge" -maxCharge "$max_precursor_charge"    

        assembly_graph_peps_mzid_psms="$working_dir"/"$metagenome_name".assembly_graph_peps.fixedKR.mzid
        java "$java_alloc" -jar "$MSGFP_jar" -s "$mgf_file" -d "$assembly_graph_peps_decoy_db" -o "$assembly_graph_peps_mzid_psms" \
        -inst "$instrument_type" -t "$parent_mass_tol" -ti "$isotope_error_range" \
        -mod "$mods_file" -tda "$is_target_decoy_db" \
        -minCharge "$min_precursor_charge" -maxCharge "$max_precursor_charge"

        # 5.b. Convert output to tsv file
        long_contig_peps_tsv_psms="$working_dir"/"$metagenome_name".long_contig_peps.fixedKR.tsv
        java "$java_alloc" -cp "$MSGFP_jar" "$MSGFP_mzid_to_tsv" \
        -i "$long_contig_peps_mzid_psms" -showDecoy "$report_decoy_psms"

        assembly_graph_peps_tsv_psms="$working_dir"/"$metagenome_name".assembly_graph_peps.fixedKR.tsv
        java "$java_alloc" -cp "$MSGFP_jar" "$MSGFP_mzid_to_tsv" \
        -i "$assembly_graph_peps_mzid_psms" -showDecoy "$report_decoy_psms"

        # 6. Retrieve results meeting FDR threshold
        long_contig_peps_fdr_psms_old="$working_dir"/"$metagenome_name".long_contig_peps.fixedKR.tsv."$fdr_level".tsv
        long_contig_peps_fdr_psms="${long_contig_peps_fdr_psms_old%.tsv."$fdr_level".tsv}"."$fdr_level".tsv
        python "$parseFDR_script" "$long_contig_peps_tsv_psms" "$fdr_level"
        mv "$long_contig_peps_fdr_psms_old" "$long_contig_peps_fdr_psms"

        assembly_graph_peps_fdr_psms_old="$working_dir"/"$metagenome_name".assembly_graph_peps.fixedKR.tsv."$fdr_level".tsv
        assembly_graph_peps_fdr_psms="${assembly_graph_peps_fdr_psms_old%.tsv."$fdr_level".tsv}"."$fdr_level".tsv
        python "$parseFDR_script" "$assembly_graph_peps_tsv_psms" "$fdr_level"
        mv "$assembly_graph_peps_fdr_psms_old" "$assembly_graph_peps_fdr_psms"

        # 7. Combine FGSP and MSGF+ PSMs into single peptide list
        combined_fdr_psms="$working_dir"/"$metagenome_name".mg_combined."$fdr_level".tsv
        python "$combineFragandDBGraph_script" "$long_contig_peps_fdr_psms" "$long_contig_pep_dna" "$long_contig_peps" \
        "$assembly_graph_peps_fdr_psms" "$assembly_graph_peps" "$long_assembly" "$combined_fdr_psms"

        # 8. Find protein sequences from assembly graph and combined PSM list
        protein_db="$working_dir"/"$metagenome_name".proteins.fasta
        "$DBGraphPep2Pro_bin" -e "$assembly_graph" -s "$assembly" -p "$combined_fdr_psms" -o "$protein_db"

        # 9.a. Spectra search against protein database
        protein_mzid_psms="$working_dir"/"$metagenome_name".proteins.mzid
        java "$java_alloc" -jar "$MSGFP_jar" -s "$mgf_file" -d "$protein_db" -o "$protein_mzid_psms" \
        -inst "$instrument_type" -t "$parent_mass_tol" -ti "$isotope_error_range" \
        -mod "$mods_file" -tda "$is_target_decoy_db" \
        -minCharge "$min_precursor_charge" -maxCharge "$max_precursor_charge"

        # 9.b. Convert output to tsv file
        protein_tsv_psms="$working_dir"/"$metagenome_name".proteins.tsv
        java "$java_alloc" -cp "$MSGFP_jar" "$MSGFP_mzid_to_tsv" \
        -i "$protein_mzid_psms" -showDecoy "$report_decoy_psms"

        # 10. Retrieve results meeting FDR threshold
        protein_fdr_psms_old="$working_dir"/"$metagenome_name".proteins.tsv."$fdr_level".tsv
        protein_fdr_psms="${protein_fdr_psms_old%.tsv."$fdr_level".tsv}"."$fdr_level".tsv
        python "$parseFDR_script" "$protein_tsv_psms" "$fdr_level"
        mv "$protein_fdr_psms_old" "$protein_fdr_psms"

    rm "$assembly" "$assembly_graph" "$long_assembly" \
    "$long_contig_peps" "$long_contig_pep_dna" "$assembly_graph_peps" \
    "$long_contig_peps_decoy_db" "$assembly_graph_peps_decoy_db" \
    "$long_contig_peps_mzid_psms" "$assembly_graph_peps_mzid_psms" \
    "$long_contig_peps_tsv_psms" "$assembly_graph_peps_tsv_psms" \
    "$long_contig_peps_fdr_psms" "$assembly_graph_peps_fdr_psms" \
    "$combined_fdr_psms" "$protein_db" "$protein_tsv_psms"