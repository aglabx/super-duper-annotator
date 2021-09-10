from os.path import join
from os.path import abspath

settings = {
    "threads": 32,
    "input_file": abspath("./test_data/NC_000913.fna"),
    "binary_folder": abspath("./prokka/binaries/linux/"),
    "work_folder": abspath("./"),
    "output_folder": abspath("./result"),
    "aragorn_genetic_code": "11",
    "barrnap_kingdom": "bac",
    "blastp_is": abspath("./prokka/db/kingdom/Bacteria/IS"),
    "blastp_amr": abspath("./prokka/db/kingdom/Bacteria/AMR"),
    "blastp_sprot": abspath("./prokka/db/kingdom/Bacteria/sprot"),
}

tools_settings = {
    "aragorn": {
        "binary": join(settings["binary_folder"], "aragorn"),
        "input": settings["input_file"],
        "output": join(settings["output_folder"], "aragorn.out"),
        "options": "-l -gc%(aragorn_genetic_code)s -w" % settings,
    },
    "barrnap": {
        "binary": join(settings["binary_folder"], "barrnap"),
        "options": "--kingdom %(barrnap_kingdom)s --threads %(threads)s --quiet" % settings,
        "output": join(settings["output_folder"], "barrnap.out"),
        "input": settings["input_file"],
    },
    "parallel": {
        "binary": "parallel",
        "input": settings["input_file"],
        "options": "--gnu --plain -j %(threads)s --block 22707 --recstart '>' --pipe" % settings,
    },
    "prodigal": {
        "binary": join(settings["binary_folder"], "prodigal"),
        "input": settings["input_file"],
        "output": join(settings["output_folder"], "prodigal_out.faa"),
        "options": "-c -m -g 11 -p single -f sco -q -a" % settings,
    },
    "blastp": {
        "binary": join(settings["binary_folder"], "blastp"),
        "is_input": settings["input_file"],
        "is_output": join(settings["output_folder"], "is.blastp"),
        "amr_output": join(settings["output_folder"], "amr.blastp"),
        "sprot_output": join(settings["output_folder"], "sprot.blastp"),
        "is_options": "-query - -db %(blastp_is)s -evalue 1e-30 -qcov_hsp_perc 90 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
        "amr_options": "-query - -db %(blastp_amr)s -evalue 9.99999999999999e-301 -qcov_hsp_perc 90 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
        "sprot_options": "-query - -db %(blastp_sprot)s -evalue 1e-09 -qcov_hsp_perc 80 \
        -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
    },
    "makeblastdb": {
        "binary": join(settings["binary_folder"], "makeblastdb"),
        "is_input": settings["blastp_is"],
        "amr_input": settings["blastp_amr"],
        "sprot_input": settings["blastp_sprot"],
        "options": "-hash_index -dbtype prot -in",
    }
}
