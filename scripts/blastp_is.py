def makeblastdb_is_command(makeblastdb_binary, makeblastdb_options, makeblastdb_is_input,):
    makeblastdb_is_command = f"{makeblastdb_binary} {makeblastdb_options} {makeblastdb_is_input}"
    return makeblastdb_is_command

#def parallel_command(parallel_binary, parallel_input, parallel_options):
    #parallel_command = f"{parallel_binary} {parallel_options}"
    #return parallel_command

#def prodigal_command(prodigal_binary, prodigal_input, prodigal_output, prodigal_options):
    #faa_input = f"{prodigal_options} {prodigal_output}"
    #return faa_input

def blastp_is_command(blastp_binary, blastp_is_input, blastp_is_output, blastp_is_options):
    blastp_is_command = f"{blastp_binary} {blastp_is_output} {blastp_is_options}"
    return blastp_is_command

def parallel_blastp_is_command(faa_input, blastp_is_command, parallel_command):
    parallel_blastp_is_command = f"cat{faa_input} | {parallel_command} {blastp_is_command}"
    return parallel_blastp_is_command
