def getInputChannels() {
	return channel.of(params.message)
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