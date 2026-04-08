FROM materialscloud/tools-barebone:1.4.0

LABEL maintainer="Osvaldo Hernandez-Cuellar <osvaldo.hernandezcuellar@epfl.ch>, Liam O. Marsh <liam.marsh@epfl.ch>, and Ruben Laplaza <ruben.laplazasolanas@epfl.ch>"

COPY ./requirements.txt /home/app/code/requirements.txt

RUN pip3 install -U pip setuptools wheel

# install NumPy 2.x first (because wheels/extensions may be selected based on it)
RUN pip3 install --no-cache-dir "numpy>=2,<3"

# FORCE replace pymatgen in the system site-packages
RUN pip3 install --no-cache-dir --upgrade --force-reinstall pymatgen

#RUN apt-get update && apt-get install -y libxrender-dev libxext-dev \
# && rm -rf /var/lib/apt/lists/*
RUN rm -f /etc/apt/sources.list.d/passenger.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends libxrender-dev libxext-dev \
 && rm -rf /var/lib/apt/lists/*

USER app
WORKDIR /home/app/code

# Install web requirements (numpy, etc.)
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

