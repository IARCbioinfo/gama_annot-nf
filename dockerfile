################## BASE IMAGE #####################
FROM continuumio/miniconda3:22.11.1

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="22.11.1"
LABEL software="gama_annot-nf"
LABEL software.version="2.0"
LABEL about.summary="Container image containing all requirements for gama_annot-nf"
LABEL about.home="http://github.com/IARCbioinfo/gama_annot-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/gama_annot-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/gama_annot-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **cahaisv** <**cahaisv@iarc.who.int**>

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n gama_annot-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/gama_annot-nf/bin:$PATH
