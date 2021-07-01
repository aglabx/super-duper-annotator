def makeblastdb_sprot_command(makeblastdb_binary, makeblastdb_options, makeblastdb_sprot_input,):
    makeblastdb_sprot_command = f"{makeblastdb_binary} {makeblastdb_options} {makeblastdb_sprot_input}"
    return makeblastdb_sprot_command

#def parallel_command(parallel_binary, parallel_input, parallel_options):
    #parallel_command = f"{parallel_binary} {parallel_options}"
    #return parallel_command

#def prodigal_command(prodigal_binary, prodigal_input, prodigal_output, prodigal_options):
    #faa_input = f"{prodigal_options} {prodigal_output}"
    #return faa_input

def blastp_sprot_command(blastp_binary, blastp_sprot_input, blastp_sprot_output, blastp_sprot_options):
    blastp_sprot_command = f"{blastp_binary} {blastp_sprot_output} {blastp_sprot_options}"
    return blastp_sprot_command

def parallel_blastp_sprot_command(faa_input, blastp_sprot_command, parallel_command):
    parallel_blastp_sprot_command = f"cat{faa_input} | {parallel_command} {blastp_sprot_command}"
    return parallel_blastp_sprot_command
