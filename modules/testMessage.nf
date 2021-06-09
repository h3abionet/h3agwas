def getInputChannels() {
    return Channel.of(params.message)
}

def getChannels() {
    return Channel.of(params.messa)
}

process printToScreen {

    input:
        val messa
        val message

    output:
        stdout
    
    script:
        """
        echo ${messa}, ${message}
        """
}