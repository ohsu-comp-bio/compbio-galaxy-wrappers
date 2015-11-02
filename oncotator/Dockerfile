FROM python:2

ENV ONCOTATOR_DIR /opt/oncotator

WORKDIR $ONCOTATOR_DIR

RUN apt-get update && apt-get install unzip

RUN pip install numpy

RUN wget --no-check-certificate 'https://github.com/elephanthunter/PyVCF/archive/master.zip' && \
	unzip master.zip && \
	cd PyVCF-master && \
	python setup.py install && \
	cd .. && \
	rm -Rf PyVCF-master && \
	rm -f master.*

RUN pip install ngslib 

RUN wget https://github.com/broadinstitute/oncotator/archive/v1.8.0.0.tar.gz && \
	tar --strip-components=1 -xvzf v1.8.0.0.tar.gz && \
 	python setup.py install

RUN easy_install numpy
