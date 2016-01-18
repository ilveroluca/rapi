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
* python 2.x
* Swig 3

Plus everything necessary to build C programa.


Running tests
------------------

Run `make tests`



Using it
---------------

That's too complicated to explain in this README file :-D

The API is defined in `include/rapi.h`.  The only implementation existing at the
momemnt wraps BWA-MEM; you can find it under `rapi_bwa`, where you'll also find
a static library (after building it, of course) that you can link to your own
programs.

The more interesting bit is probably the Python interface to RAPI, which you can
find under `pyrapi`.  If you have all the necessary


FAQ
-------



### Is the API stable?

No, it's not.  The implementation works, but the API is still under development
and thus will change.  If you use it, do provide your suggestions on how to
improve it.
