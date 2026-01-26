FROM materialscloud/tools-barebone:1.4.0

LABEL maintainer="Osvaldo Hernandez-Cuellar <osvaldo.hernandezcuellar@epfl.ch>, Liam O. Marsh <liam.marsh@epfl.ch>, and Ruben Laplaza <ruben.laplazasolanas@epfl.ch>"

COPY ./requirements.txt /home/app/code/requirements.txt

RUN pip3 install -U 'pip>=10' setuptools==65.4.1 wheel

RUN apt-get update && apt-get install -y libxrender-dev libxext-dev \
 && rm -rf /var/lib/apt/lists/*

USER app
WORKDIR /home/app/code

# Install web requirements (numpy==1.26.4 etc.)
RUN pip3 install --user -r requirements.txt

USER root

COPY ./config.yaml /home/app/code/webservice/static/config.yaml
COPY ./user_templates/ /home/app/code/webservice/templates/user_templates/
COPY ./user_static/ /home/app/code/webservice/user_static/
COPY ./compute/ /home/app/code/webservice/compute/
COPY ./web_module.py /home/app/code/webservice/
COPY ./base_templates/* /home/app/code/webservice/templates/
COPY ./toolsbarebone_mod/__init__.py /home/app/code/tools_barebone/structure_importers/

RUN chown -R app:app /home/app/code/webservice/

