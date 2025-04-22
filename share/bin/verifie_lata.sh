#!/bin/bash

# Script for lata comparison used in Front tracking test cases
# Done because interface cannot be written in lml format

[ ! -f $1.lata ] && echo "File $1.lata is missing. Test Failed" && exit 1


# Trouve le nom du cas test:
name=$1

echo "Test case name: ${name}"

# untar ref for compare_lata
tar xzf ref.tar.gz


# Update ref if asked
if [ "$UPDATE_LATA_REF" != "" ]; then

	test_path=${TrioCFD_project_directory}/tests/Reference/$(realpath -s --relative-to=$TRUST_TMP/tests $(pwd))

    echo $(pwd)
    echo ${test_path}
    
    echo "Updating ref"
    rm -r ref
    mkdir -p ref
    cp ${name}.lata* ref/
    tar czf ref.tar.gz ref

    # check ref size (max 50 MB)
	max_size=50000000
	ref_size=$(wc -c <"ref.tar.gz")
	if [ $ref_size -ge $max_size ]; then
		echo "Error: reference size too big : $ref_size"
		exit 1
	fi

    cp ref.tar.gz ${test_path}
fi


# Compare to ref
for cas in "" "PAR_"
do
	if [ -f ${cas}${name}.lata ] ; then 
		echo Comparing ${cas}${name}.lata ref/${name}.lata
		if compare_lata ${cas}${name}.lata ref/${name}.lata 
		then
			echo "ok"
		else
			echo "Comparaison avec le lata de reference echouee."
			exit 1
		fi
	fi
done
