#!/bin/bash

#exit if any command fails
set -e

if [ `whoami` != "root" ] ; then
    echo "Rerun as root"
    exit 1
fi

apt -y update
apt -y upgrade
apt-get install -y fftw3-dev libhdf5-dev libhdf5-serial-dev python3-pip build-essential

