process printMessage {
	input:
		val message
	output:
		stdout
	script:
		"""
		echo ${message}
		"""
}