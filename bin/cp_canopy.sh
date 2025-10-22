#!/bin/bash

for proc in $(ls -d proc*); do
	for t in $(ls -d processor0/[1-9]* | cut -d "/" -f 2); do
		cp -f $proc/0/canopy_* $proc/$t/ 
	done
	echo $proc
done
