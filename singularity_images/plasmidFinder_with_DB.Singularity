Bootstrap: docker
From: ubuntu:16.04

%files

%environment

   PATH=$PATH:/opt/plasmidfinder:/blast+/ncbi-blast-2.10.0+/bin

%labels

   AUTHOR nvlachos

%files
 /scicomp/home/nvx4/Singularity/ncbi-blast-2.10.0+-x64-linux.tar.gz /

%post
   #install ubuntu dependencies
   apt-get update
   apt-get install -y build-essential python3 git libz-dev python3-pip

   #install plasmidFinder
   cd /opt
   git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
   git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
   git clone https://bitbucket.org/genomicepidemiology/kma.git
   cd kma && make
   cd /opt/plasmidfinder_db
   python3 INSTALL.py ../kma/kma_index

   pip3 install --upgrade pip
   pip3 install tabulate
   pip3 install cgecore
   pip3 install Biopython

   #install blast+
   mkdir /blast+
   mv /ncbi-blast-2.10.0+-x64-linux.tar.gz /blast+
   chmod -R 755 /blast+
   cd /blast+
   tar zxvpf ncbi-blast-2.10.0+-x64-linux.tar.gz

   #clean-up
   apt-get clean

%runscript

