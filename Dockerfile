FROM materialscloud/tools-barebone:1.3.0

LABEL maintainer="Developer Name <developer.email@xxx.yy>"

# Python requirements
COPY ./requirements.txt /home/app/code/requirements.txt
COPY ./cell2mol/requirements.txt /home/app/code/compute-requirements.txt
# Run this as sudo to replace the version of pip

RUN pip3 install -U 'pip>=10' setuptools wheel
# install packages as normal user (app, provided by passenger)

RUN apt-get update
RUN apt-get install -y libxrender-dev libxext-dev

USER app
WORKDIR /home/app/code
# Install pinned versions of packages
COPY ./requirements.txt /home/app/code/requirements.txt
COPY ./cell2mol/requirements.txt /home/app/code/compute-requirements.txt
RUN pip3 install --user -r requirements.txt
RUN pip3 install --user -r compute-requirements.txt

# Go back to root.
# Also, it should remain as user root for startup
USER root

# Copy various files: configuration, user templates, the actual python code, ...
COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./user_static/ /home/app/code/webservice/user_static/
COPY ./cell2mol/ /home/app/code/webservice/cell2mol/
COPY ./compute/ /home/app/code/webservice/compute/

# If you put any static file (CSS, JS, images),
#create this folder and put them here
# COPY ./user_static/ /home/app/code/webservice/user_static/

###
# Copy any additional files needed into /home/app/code/webservice/
###

# Set proper permissions on files just copied
RUN chown -R app:app /home/app/code/webservice/

USER app
WORKDIR /home/app/code/webservice/cell2mol/
RUN python3 ./setup.py build
USER root
