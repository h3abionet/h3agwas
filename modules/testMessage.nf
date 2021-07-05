def getInputChannels() {
	return channel.of(params.testMessage)
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