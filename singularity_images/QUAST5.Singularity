Bootstrap: docker 
From: ubuntu:16.04

%labels 
	maintainer Yoann Dufresne 
	<yoann.dufresne@pasteur.fr> quast.version 5.0.0
%post
	##### System #####
	apt update -y 
	apt upgrade -y 
	apt install -y build-essential wget zlib1g-dev pkg-config libfreetype6-dev libpng-dev python-matplotlib perl openjdk-8-jdk python3 python-setuptools python-dev language-pack-en-base
	
	##### Quast #####
	wget -P /tmp/ https://downloads.sourceforge.net/project/quast/quast-5.0.0.tar.gz 
	cd /tmp 
	tar xvzf /tmp/quast-5.0.0.tar.gz 
	rm /tmp/quast-5.0.0.tar.gz 
	if [ -d "/quast" ]; then
		rm -rf /quast/ 
	fi 
	mv /tmp/quast-5.0.0 /quast 
	cd /quast/ 
	./install.sh

%runscript 
	python3 /quast/quast.py $@
