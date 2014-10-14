Full information regarding the Cosmic clocks data can be found here

http://www.physics-astronomy.unibo.it/en/research/areas/astrophysics/cosmology-with-cosmic-chr
onometers

Three files are provided, but note that only one of them must be used
for any given run (i.e they are NOT independent data):
- the "BC03" file report the H(z) measurements from Moresco et al.
  (2012) analysis based on Bruzual & Charlot (2003) stellar population
  synthesis models
- the "MaStro" file report the H(z) measurements from Moresco et al.
  (2012) analysis based on Maraston & Stromback (2011) stellar
  population synthesis models
- the "BC03_all" file report the H(z) measurements contained in "BC03"
  file combined with [Simon et al. 2005, PhRvD, 71, 123001; Stern et
  al. 2010, JCAP, 2, 8] measurements.

In terms of Monte Python, this means that you can specify either one
of the three new experiments in your experiment list
- 'cosmic_clocks_BC03'
- 'cosmic_clocks_BC03_all'
- 'cosmic_clocks_MaStro'
If you specify two or more, an error will be raised.

The relevant publications are:
- for the data at z<1.1 (i.e files with “BC03” or “MaStro” in the
  name) Moresco et al. 2012 JCAP , 08, 006
- for the data up to z\sim1.7  (i.e. the files with “BC03_all” in the
  name) then also look at Simon et al. 2005, PhRvD, 71, 123001; Stern
  et al. 2010, JCAP, 2, 8

Since all the other previous analysis in literature (Simon et al.
2005, Stern et al. 2010) provided H(z) measurements only for BC03
models, for the sake of uniformity, in the last file are combined
estimates obtained with BC03 models. While the dependence on stellar
population synthesis model is unimportant for most cosmological
analyses (as shown in Moresco et al.  2012), and "BC03_all" can be
safely used, if you want to be conservative the “BC03” and “MaStro”
files are also provided.
