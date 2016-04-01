rapi
====

Read aligner API

This project defines an API for a short read aligner (to map short reads to a
reference sequence).  It also includes a reference implementation that wraps
BWA-MEM and Python bindings for the API.


Compiling
---------------

Run `make`.

You'll need to have some tools on your system:

* git
* Swig 3
* python 2.x, including dev files
* python-setuptools
* zlib-dev

Plus everything necessary to build C programs.


Running tests
------------------

Run `make tests`



Using it
---------------

That's too complicated to explain in this README file :-D

C interface
+++++++++++++++

The API is defined in `include/rapi.h`.  You can find an example showing how to
use it under `example/rapi_test.c`.  The only plug-in implementation existing at
the moment wraps BWA-MEM; you can find it under `rapi_bwa`, where you'll also
find a static library (after building it, of course) that you can link to your
own programs.


Python interface
+++++++++++++++++++

An equally interesting bit is the Python interface to RAPI, which you can
find under `pyrapi`.  The script `example/align.py` implements a command-line
interface to the aligner (you can use it to align reads in fastq files and
generate SAM).


Authors
---------

RAPI is written by Luca Pireddu (pireddu@crs4.it).

Thanks go to these people for their contributions:
  * Riccardo Berutti (berutti@berutti.net)
  * Sebastian Schoenherr (sebastian.schoenherr@uibk.ac.at)


FAQ
-------



### Is the API stable?

No, it's not.  The implementation works, but the API is still under development
and thus will change.  If you use it, do provide your suggestions on how to
improve it.
