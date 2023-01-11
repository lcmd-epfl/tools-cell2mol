##############
tools-cell2mol
##############

In this repository we provide the code to deploy an online service for
the interpretation of crystallographic unit cells using `cell2mol_`.
A `live demo`_ is currently hosted on the `Materials Cloud`_ web portal.

============
Contributors
============

Developed by LCMD members at EPFL under the supervision of Prof. Corminboeuf:

- Sergi Vela, Yuri Cho, Ruben Laplaza (EPFL) [cell2mol development]
- Liam Marsh, Osvaldo Cuellar (EPFL) [tools-cell2mol development]

=========================================
Docker image and running the tool locally
=========================================
Docker images are automatically built and hosted on `DockerHub under the repository materialscloud/tools-cell2mol`_.

If you want to run locally the latest version, you can execute::

  docker pull materialscloud/tools-cell2mol:latest
  docker run -p 8093:80 materialscloud/tools-cell2mol:latest

and then connect to ``http://localhost:8090`` with your browser.


.. _Materials Cloud: https://www.materialscloud.org/work/tools/cell2mol
.. _cell2mol: https://github.com/lcmd-epfl/cell2mol
.. _DockerHub under the repository materialscloud/tools-cell2mol: https://hub.docker.com/repository/docker/materialscloud/tools-cell2mol
.. _live demo: https://cell2mol.materialscloud.io/

