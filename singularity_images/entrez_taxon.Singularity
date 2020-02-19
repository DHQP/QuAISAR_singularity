
Bootstrap: docker
From: ubuntu:16.04

%files

%environment

   #Add gottcha bin folder
   PATH=$PATH:/GAMA:/usr/bin:/bin

%labels

   AUTHOR nvlachos

%files
   /scicomp/home/nvx4/Singularity/entrez_get_taxon_from_accession.py /
   
%post
   #install ubuntu dependencies
   apt-get update   
   apt-get  -y upgrade
   apt-get install -y build-essential python3 python3-pip python-pip wget unzip libz-dev libncurses5-dev libkrb5-3

   #install python dependencies
   pip3 install --upgrade pip
   pip3 install biopython

   mkdir /entrez
   mv /entrez_get_taxon_from_accession.py /entrez
   chmod -R 755 /entrez

   #clean up
   apt-get clean
   rm -rf /var/lib/apt/lists/*

%runscript

