#!/bin/bash

linkFrom=$1
name=$2
linkTo=$3

if [[ ! -d "${linkTo}/${name}" ]]; then
    ln -s ${linkFrom} ${linkTo}/${name}
    echo ${linkFrom}
    echo ${linkTo}/${name}
fi
