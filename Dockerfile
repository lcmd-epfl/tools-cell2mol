FROM materialscloud/tools-barebone:1.3.0

LABEL maintainer="Osvaldo Hernandez-Cuellar <osvaldo.hernandezcuellar@epfl.ch>, Liam O. Marsh <liam.marsh@epfl.ch>, and Ruben Laplaza <ruben.laplazasolanas@epfl.ch>"

# Python requirements
COPY ./requirements.txt /home/app/code/requirements.txt
# Run this as sudo to replace the version of pip

RUN pip3 install -U 'pip>=10' setuptools wheel
# install packages as normal user (app, provided by passenger)

RUN apt-get update
RUN apt-get install -y libxrender-dev libxext-dev

USER app
WORKDIR /home/app/code
# Install pinned versions of packages
COPY ./requirements.txt /home/app/code/requirements.txt
RUN pip3 install --user -r requirements.txt

# Go back to root.
# Also, it should remain as user root for startup
USER root

# Copy various files: configuration, user templates, the actual python code, ...
COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./user_static/ /home/app/code/webservice/user_static/
COPY ./compute/ /home/app/code/webservice/compute/
#Needed to allow only .cif file formats in the upload_structure_block
COPY ./web_module.py /home/app/code/webservice/ 
COPY ./base_templates/* /home/app/code/webservice/templates/

# If you put any static file (CSS, JS, images),
#create this folder and put them here

###
# Copy any additional files needed into /home/app/code/webservice/
###

# Set proper permissions on files just copied
RUN chown -R app:app /home/app/code/webservice/


