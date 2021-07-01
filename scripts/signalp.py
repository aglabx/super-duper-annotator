def signalp_command(signalp_binary, signalp_input, signalp_output, signalp_options):
    signalp_command = f"{signalp_options} {signalp_output} 2> /dev/null"
    return signalp_command
