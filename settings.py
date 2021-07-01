import os
import os.path

settings = {
    "threads": 32,
    "input_file": os.path.abspath("./test_data/NC_000913.fna"),
    "binary_folder": os.path.abspath("./prokka/binaries/linux/"),
    "work_folder": os.path.abspath("./"),
    "output_folder": os.path.abspath("./result"),
    "aragorn_genetic_code": "11",
    "barrnap_kingdom": "bac",
    "blastp_is" : os.path.abspath("./prokka/db/kingdom/Bacteria/IS"),
    "blastp_amr" : os.path.abspath("./prokka/db/kingdom/Bacteria/AMR"),
    "blastp_sprot" : os.path.abspath("./prokka/db/kingdom/Bacteria/sprot"),
}

tools_settings = {
    "aragorn_options": {
        "aragorn_binary": os.path.join(settings["binary_folder"], "aragorn"),
        "aragorn_input": settings["input_file"],
        "aragorn_output": os.path.join(settings["output_folder"], "aragorn.out"),
        "aragorn_options": "-l -gc%(aragorn_genetic_code)s -w" % settings,
    },
    "barrnap_options": {
        "barrnap_binary": os.path.join(settings["binary_folder"], "barrnap"),
        "barrnap_options": "--kingdom %(barrnap_kingdom)s --threads %(threads)s --quiet" % settings,
        "barrnap_output": os.path.join(settings["output_folder"], "barrnap.out"),
        "barrnap_input": settings["input_file"],
    },
    "parallel_options": {
        "parallel_binary": "parallel",
        "parallel_input": settings["input_file"],
        "parallel_options": "--gnu --plain -j %(threads)s --block 22707 --recstart '>' --pipe" % settings,
    },
    "prodigal_options": {
        "prodigal_binary": os.path.join(settings["binary_folder"], "prodigal"),
        "prodigal_input": settings["input_file"],
        "prodigal_output": os.path.join(settings["output_folder"], "prodigal_out.faa"),
        "prodigal_options": "-c -m -g 11 -p single -f sco -q -a" % settings,
    },
    "blastp_options": {
        "blastp_binary": os.path.join(settings["binary_folder"], "blastp"),
        "blastp_is_input": settings["input_file"],
        "blastp_is_output": os.path.join(settings["output_folder"], "is.blastp"),
        "blastp_amr_output": os.path.join(settings["output_folder"], "amr.blastp"),
        "blastp_sprot_output": os.path.join(settings["output_folder"], "sprot.blastp"),
        "blastp_is_options": "-query - -db %(blastp_is)s -evalue 1e-30 -qcov_hsp_perc 90 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
        "blastp_amr_options": "-query - -db %(blastp_amr)s -evalue 9.99999999999999e-301 -qcov_hsp_perc 90 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
        "blastp_sprot_options": "-query - -db %(blastp_sprot)s -evalue 1e-09 -qcov_hsp_perc 80 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
    },
    "makeblastdb_options": {
        "makeblastdb_binary": os.path.join(settings["binary_folder"], "makeblastdb"),
        "makeblastdb_is_input": settings["blastp_is"],
        "makeblastdb_amr_input": settings["blastp_amr"],
        "makeblastdb_sprot_input": settings["blastp_sprot"],
        "makeblastdb_options": "-hash_index -dbtype prot -in",    
    },
     "signalp_options": {
        "signalp_binary": os.path.join(settings["binary_folder"], "signalp"),
        "signalp_input": settings["input_file"],
        "signalp_output": os.path.join(settings["output_folder"], "signalp.out"),
        "signalp_options": "-tmp ./result -prefix ./result/signalp -org gram+ -format short -fasta" % settings,
    }
}