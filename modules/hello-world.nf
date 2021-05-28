process printMessage {
    tag "${message}"

	input:
		val message
	output:
		stdout
	script:
		"""
		echo ${message}
		"""
}
