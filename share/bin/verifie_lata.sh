#!/bin/bash

# Script for lata comparison used in Front tracking test cases
# Done because interface cannot be written in lml format

[ ! -f $1.lata ] && echo "File $1.lata is missing. Test Failed" && exit 1


# Trouve le nom du cas test:
fullname=$1

name=${fullname#PAR_}

echo "Test case name: ${fullname}"

# untar ref for compare_lata
tar xzf ref.tar.gz


# Update ref if asked
if [ "$UPDATE_LATA_REF" != "" ] && [ "$fullname" != "$name" ]
then

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
if [ -f ${cas}${name}.lata ] ; then 

    # differences between seq and parallel are allowed if # ECART_SEQ_PAR # is found in the datafile
    ecart_seq_par=$(grep -n "ECART_SEQ_PAR" ${name}.data | tail -1)
    ecart_seq_par=${ecart_seq_par%%:*}
    if [ "$fullname" != "$name" ] && [ "$ecart_seq_par" != "" ]
    then
        echo "Warning: allowing SEQ/PAR differences according to ${name}.data, at line ${ecart_seq_par}."
    fi

    # comparison
    echo Comparing ${fullname}.lata ref/${name}.lata
    if compare_lata ${fullname}.lata ref/${name}.lata 
    then
        echo "ok"
        exit 0
    else

        if [ "$fullname" != "$name" ] && [ "$ecart_seq_par" != "" ]
        then
            echo "Warning: found SEQ/PAR differences, but this was expected."
            exit 0
        else
        echo "Error: Failed comparison to reference lata."
        exit 1
        fi
        
    fi
fi
