Bootstrap: docker
From: ubuntu:16.04

%files

%environment

   PATH=$PATH:/opt/bowtie2-2.1.0:/opt/samtools-0.1.18

%labels

   AUTHOR jyao


%post
   #install ubuntu dependencies
   apt-get update
   apt-get install -y build-essential python-pip git wget unzip libz-dev libncurses5-dev

   #install python dependency
   pip install --upgrade pip
   pip install scipy

   #install samtools
   cd /opt
   wget https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
   tar xjf samtools-0.1.18.tar.bz2
   cd samtools-0.1.18
   make
#   make install
   cd ..
#   apt-cache madison samtools
#   apt-get install samtools=0.1.18-1ubuntu2


   #install bowtie2
   wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip
   unzip bowtie2-2.1.0-linux-x86_64.zip 

   #install SRST2
   pip install git+https://github.com/katholt/srst2

   #clean up
   rm /opt/samtools-0.1.18.tar.bz2
   rm /opt/bowtie2-2.1.0-linux-x86_64.zip
   apt-get clean
   rm -rf /var/lib/apt/lists/*

%runscript

