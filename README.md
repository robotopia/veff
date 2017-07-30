# README #

### Overview ###

* veff is a collection of programs used primarily for the analysis of pulsar secondary spectra.
* Version: 0.8.0
* Includes the programs veff, parabfit, and a template geogebra XML file for producing visualisations of the Earth/pulsar systems

### Installation ###

* If you intend to use the veff program, make sure you have first installed the [CSPICE Tooklit for C](https://naif.jpl.nasa.gov/naif/toolkit_C.html) (not required for parabfit).
* Open the Makefile in the base directory and edit the CSPICE\_DIR and INSTALL\_DIR variables as needed.
* Run **make** and **make install**
* For more comprehensive documentation, run **make documentation** (documentation will appear in the doc/ subdirectory as doc.pdf). This requires pdflatex to be installed on your system.