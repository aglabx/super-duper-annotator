def makeblastdb_command(makeblastdb_binary, makeblastdb_options, makeblastdb_input,):
    makeblastdb_command = f"{makeblastdb_binary} {makeblastdb_options} {makeblastdb_input}"
    return makeblastdb_command