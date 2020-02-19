Bootstrap: docker
From: ubuntu:12.04

%files

%environment

   #Add gottcha bin folder
   PATH=$PATH:/opt/gottcha/bin/

%labels

   AUTHOR jyao


%post
   #install ubuntu dependencies
   apt-get update
   apt-get install -y perl build-essential git wget unzip libz-dev

   #install gottcha
   cd /opt
   git clone https://github.com/LANL-Bioinformatics/GOTTCHA.git gottcha
   cd gottcha
   ./INSTALL.sh

   #clean up
   apt-get clean
   rm -rf /var/lib/apt/lists/*

%runscript

