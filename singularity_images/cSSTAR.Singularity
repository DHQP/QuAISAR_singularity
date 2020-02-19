
Bootstrap: docker
From: ubuntu:16.04

%files

%environment

   #Add gottcha bin folder
   PATH=$PATH:/cSSTAR:/usr/bin:/blast+/ncbi-blast-2.10.0+/bin

%labels

   AUTHOR nvlachos

%files
   /scicomp/home/nvx4/Singularity/c-SSTAR_gapped.py /
   /scicomp/home/nvx4/Singularity/c-SSTAR_ungapped.py /
   /scicomp/home/nvx4/Singularity/ncbi-blast-2.10.0+-x64-linux.tar.gz /

%post
   #install ubuntu dependencies
   apt-get update   
   apt-get  -y upgrade
   apt-get install -y build-essential python3 python3-pip python-pip wget unzip libz-dev libncurses5-dev libkrb5-3

   #install python dependencies
   pip3 install --upgrade pip
   pip3 install biopython

   mkdir /cSSTAR
   mkdir /blast+
   mv /c-SSTAR_ungapped.py /cSSTAR
   mv /c-SSTAR_gapped.py /cSSTAR
   mv /ncbi-blast-2.10.0+-x64-linux.tar.gz /blast+
   chmod -R 755 /cSSTAR
   chmod -R 755 /blast+

   cd /blast+
   tar zxvpf ncbi-blast-2.10.0+-x64-linux.tar.gz

   #clean up
   apt-get clean
   rm -rf /var/lib/apt/lists/*

%runscript

