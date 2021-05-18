#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 *	PRINT HELLO WORLD
 *	=================
 *
 *	This is a test pipeline. Use it to test that you have set up your
 *	project space correctly.
 *
 *********************************************************************/

include {
	printMessage;
} from "${params.modulesDir}/hello-world.nf"

workflow {

	string = channel.from("Hello", "World!")

	printMessage(string) | view()

}
