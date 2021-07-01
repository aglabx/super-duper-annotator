def prodigal_command(prodigal_binary, prodigal_input, prodigal_output, prodigal_options):
    prodigal_command = f"{prodigal_binary} {prodigal_input} {prodigal_options} {prodigal_output}"
    faa_input = prodigal_output
    return faa_input