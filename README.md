# A contribution to the Data Challenge 1
### by Pragati Mitra and Lech Wiktor Piotrowski

Call efield2adc_tracebytrace.py to create voltage.root from Coarse2.root in this directory.

It assumes that:
* There is Coarse2.root (a file with traces from Data Challenge 1) in this directory (or you can give a path to it as a first parameter)
* There are files necessary by XDU electronic chain in the electronic_chain/XDU_electronic_chain/XDU_files directory. The files and directories should include:
  * antennaVSWR
  * cableparameter
  * filterparameter
  * LNASparameter
  * 30_250galactic.mat
  * Complex_RE.mat
