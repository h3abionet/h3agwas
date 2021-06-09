def getInputChannels() {
    return Channel.of(params.message)
}

process printToScreen {

    input:
        val message

    output:
        stdout
    
    script:
        """
        echo ${message}
        """
}