#!/usr/bin/env bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b

# Didnt work on my account CERES
echo 'export PATH="~/miniconda3/bin:$PATH"' >> ~/.bashrc

sleep 2

. ~/.bashrc
conda init bash
conda config --set auto_activate_base true
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c conda-forge singularity
rm Miniconda3-latest-Linux-x86_64.sh
