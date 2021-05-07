#!/usr/bin/env bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b
echo 'export PATH="$1/miniconda3/bin:$PATH"' >> ~/.bashrc
. ~/.bashrc
conda init bash
conda config --set auto_activate_base false
. ~/.bashrc
rm Miniconda3-latest-Linux-x86_64.sh
