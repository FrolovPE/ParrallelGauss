#!/bin/bash

if [ $# != 3 ]; then
	echo "Usage: supertest.sh <a.out> <r> <n(до какого n делать тесты)>"
else
    exe=$1
    r=$2
    n=$3
    
    if [ -f $exe ];then
    echo
    else
        echo "$exe doesn't exist"
        exit
    fi

<<<<<<< HEAD
    for ((i=2; i<=$n; i++)); do

        for ((j=1; j<i; j++)); do

            for ((k=1; k<4;k++)); do

		for ((q = 1; q <=4 ; q++)); do
                    echo "--------------------------------------------------------------------------"
                    echo "$exe $i $j $q $r $k"
                    $exe $i $j $q $r $k
                    echo ""
		done

=======
    for((i=2; i<=$n; i++)); do

        for((j=1; j<i; j++)); do

            for((k=1; k<4;k++)); do
                    echo "--------------------------------------------------------------------------"
                    echo "$exe $i $j $r $k"
                    $exe $i $j $r $k
                    echo ""
>>>>>>> e6b7ae0 ('Thu Dec 11 23:21:30 ')
            done

        done

    done

<<<<<<< HEAD
fi
=======
fi
>>>>>>> e6b7ae0 ('Thu Dec 11 23:21:30 ')
