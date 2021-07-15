FROM ribocleaner-core:latest

USER root

# Install python
RUN apt-get update && apt-get install -y python3 python3-dev python3-pip imagemagick graphviz libgraphviz-dev && \
        ln -s /usr/bin/python3 /opt/bin/python && ln -s /usr/bin/pip3 /opt/bin/pip

# Install required python packages
# pin to past version of snakemake because current version requires a special version of imagemagick
RUN pip install pandas snakemake==6.4.1 matplotlib jinja2 networkx pygments && \
	pip install pygraphviz

# make a base analysis directory and directories to store the results
# analysis, workflow, and config will be mounted from the host system
RUN mkdir /analysis
WORKDIR /analysis
RUN mkdir inputs results workflow

# add the snakemake workflow
# this runs from the context of the workflow directory
ADD Snakefile ./workflow/
ADD schemas ./workflow/schemas
ADD rules ./workflow/rules
ADD report ./workflow/report

# make the cleaner user the owner of all associated files
RUN chown -R ribocleaner:ribocleaner /analysis

USER ribocleaner
