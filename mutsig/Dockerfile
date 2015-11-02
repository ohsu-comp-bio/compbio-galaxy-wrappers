FROM phusion/baseimage

ENV MUTSIG_DIR /opt/mutsig
ENV MCR_DIR /opt/mcr

WORKDIR $MCR_DIR

RUN apt-get -qq update

RUN apt-get install -y \
	wget \
	unzip \
	xorg \
	openbox \
	xauth

RUN wget http://www.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip && unzip MCR_R2013a_glnxa64_installer.zip && ./install -mode silent -agreeToLicense yes

WORKDIR $MUTSIG_DIR

RUN wget https://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSigCV_1.4.zip && unzip MutSigCV_1.4.zip 
