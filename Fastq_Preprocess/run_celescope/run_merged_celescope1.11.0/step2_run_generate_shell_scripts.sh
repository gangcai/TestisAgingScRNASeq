#!/bin/bash
for sample in TestisMix
do
	echo $sample
	conda run -n celescope1.11.0 bash ./generate_shell_scripts.sh ${sample}
done
