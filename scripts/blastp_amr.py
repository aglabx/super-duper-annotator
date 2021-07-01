def makeblastdb_amr_command(makeblastdb_binary, makeblastdb_options, makeblastdb_amr_input,):
    makeblastdb_amr_command = f"{makeblastdb_binary} {makeblastdb_options} {makeblastdb_amr_input}"
    return makeblastdb_amr_command

#def parallel_command(parallel_binary, parallel_input, parallel_options):
    #parallel_command = f"{parallel_binary} {parallel_options}"
    #return parallel_command

#def prodigal_command(prodigal_binary, prodigal_input, prodigal_output, prodigal_options):
    #faa_input = f"{prodigal_options} {prodigal_output}"
    #return faa_input

def blastp_amr_command(blastp_binary, blastp_amr_input, blastp_amr_output, blastp_amr_options):
    blastp_amr_command = f"{blastp_binary} {blastp_amr_output} {blastp_amr_options}"
    return blastp_amr_command

def parallel_blastp_amr_command(faa_input, blastp_amr_command, parallel_command):
    parallel_blastp_amr_command = f"cat{faa_input} | {parallel_command} {blastp_amr_command}"
    return parallel_blastp_amr_command
