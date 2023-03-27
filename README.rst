##############
tools-cell2mol
##############

In this repository we provide the code to deploy an online service for
the interpretation of crystallographic unit cells using `cell2mol_`.

.. A `live demo`_ is currently hosted on the `Materials Cloud`_ web portal.

.. contents::

.. section-numbering::   

===========
How to cite
===========

Vela, S.; Laplaza, R.; Cho, Y.; Corminboeuf, C. cell2mol: encoding chemistry to interpret crystallographic data. *Npj Comput. Mater.*
**2022**. *8*. 188; `doi: 10.1038/s41524-022-00874-9 <https://doi.org/10.1038/s41524-022-00874-9>`_


============
Contributors
============

Developed by LCMD members at EPFL under the supervision of Prof. Corminboeuf:

- Sergi Vela, Yuri Cho, Ruben Laplaza (EPFL) [cell2mol development]
- Liam Marsh, Osvaldo Hernandez-Cuellar (EPFL) [tools-cell2mol development]


=============
Prerequisites
=============

- `docker`_ >= v18.09
- `tools-barebone`_  >= 1.0.0

========================
Running the tool locally
========================

Download the `tools-cell2mol`_ repository.

In the main directory, build the docker container by executing::

    sudo docker-compose build

this has to be done only the first time or when changing the `tools-cell2mol`_ version.


To run `tools-cell2mol`_ locally, you can execute::

    sudo docker-compose up

in the main directory, and then connect to ``http://localhost:8090`` (or the direction docker indicates you in case you have 
more than one container running) with your browser.


============
Known errors
============

Few chemical patterns tend to be poorly interpreted by cell2mol, because of inconsistencies in either cell2mol itself, or in xyz2mol.

cell2mol determines the bond order between atoms based on their connectivity, typical atomic valence electrons of the atoms 
involved, and the most plausible total charge of the molecule. cell2mol is inevitably incorrect if there is an extra electron 
or radical chemical species. Other known ligands with wrong interpretations are the triiodide (I3-), and azide (N3-) ions. 
Future development of cell2mol will aim a fixing those errors. If users identify other common misinterpretations, please contact 
the authors.


.. Not yet
.. =========================================
.. Docker image and running the tool locally
.. =========================================
.. Docker images are automatically built and hosted on `DockerHub under the repository materialscloud/tools-cell2mol`_.
.. 
.. If you want to run locally the latest version, you can execute::
.. 
..   docker pull materialscloud/tools-cell2mol:latest
..   docker run -p 8093:80 materialscloud/tools-cell2mol:latest
.. 
.. and then connect to ``http://localhost:8090`` with your browser.


.. _docker: https://www.docker.com/
.. _tools-barebone: https://github.com/materialscloud-org/tools-barebone
.. _tools-cell2mol: https://github.com/lcmd-epfl/tools-cell2mol
.. _Materials Cloud: https://www.materialscloud.org/work/tools/cell2mol
.. _cell2mol: https://github.com/lcmd-epfl/cell2mol
.. _DockerHub under the repository materialscloud/tools-cell2mol: https://hub.docker.com/repository/docker/materialscloud/tools-cell2mol
.. _live demo: https://cell2mol.materialscloud.io/

