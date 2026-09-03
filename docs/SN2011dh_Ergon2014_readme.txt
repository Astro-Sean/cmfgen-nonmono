J/A+A/562/A17       SN 2011dh - The first 100 days               (Ergon+, 2014)
================================================================================
Optical and near-infrared observations of SN 2011dh - The first 100 days.
    Ergon M., Sollerman J., Fraser M., Pastorello A., Taubenberger S.,
    Elias-Rosa N., Bersten M., Jerkstrand A., Benetti S., Botticella M.T.,
    Fransson C., Harutyunyan A., Kotak R., Smartt S., Valenti S., Bufano F.,
    Cappellaro E., Fiaschi M., Howell A., Kankare E., Magill L., Mattila S.,
    Maund J., Naves R., Ochner P., Ruiz J., Smith K., Tomasella L., Turatto M.
   <Astron. Astrophys. 562, A17 (2014)>
   =2014A&A...562A..17E
================================================================================
ADC_Keywords: Supernovae ; Photometry, UBVRI ; Photometry, ugriz ;
              Photometry, infrared
Keywords: supernovae: general - supernovae: individual: SN 2011dh -
          galaxies: clusters: individual: M51

Abstract:
    We present optical and near-infrared (NIR) photometry and spectroscopy
    of the Type IIb supernova (SN) 2011dh for the first 100 days. We
    complement our extensive dataset with Swift ultra-violet (UV) and
    Spitzer mid-infrared (MIR) data to build a UV to MIR bolometric
    lightcurve using both photometric and spectroscopic data.
    Hydrodynamical modelling of the SN based on this bolometric lightcurve
    have been presented in Bersten et al. (2012ApJ...757...31B).

Description:
    Johnson-Cousins UBVRI, SDSS ugriz and 2MASS JHK magnitudes and
    pseudo-bolometric UV-MIR bolometric luminosity for SN 2011dh during
    the first 100 days.

Objects:
    ---------------------------------------------------
       RA   (2000)   DE        Designation(s)
    ---------------------------------------------------
    13 30 05.12  +47 10 11.55  SN 2011dh = SN 2011dh
    ---------------------------------------------------

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table3.dat    114       90   Johnson-Cousins UBVRI magnitudes for SN 2011dh
table5.dat    114       45   SDSS ugriz magnitudes for SN 2011dh
table6.dat     82       23   2MASS JHK magnitudes for SN 2011dh
table10.dat    48       97   Bolometric luminosity for SN 2011dh
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8  F8.2  d       Epoch     Epoch (JD-2400000)
  11- 16  F6.2  d       Phase     Days since explosion
  19- 24  F6.3  mag     Umag      ? Johnson-Cousins U magnitude (1)
  27- 32  F6.3  mag   e_Umag      ? Johnson-Cousins U magnitude error (1)
  35- 40  F6.3  mag     Bmag      ? Johnson-Cousins B magnitude (1)
  43- 48  F6.3  mag   e_Bmag      ? Johnson-Cousins B magnitude error (1)
  51- 56  F6.3  mag     Vmag      ? Johnson-Cousins V magnitude (1)
  59- 64  F6.3  mag   e_Vmag      ? Johnson-Cousins V magnitude error (1)
  67- 72  F6.3  mag     Rmag      ? Johnson-Cousins R magnitude (1)
  75- 80  F6.3  mag   e_Rmag      ? Johnson-Cousins R magnitude error (1)
  83- 88  F6.3  mag     Imag      ? Johnson-Cousins I magnitude (1)
  91- 96  F6.3  mag   e_Imag      ? Johnson-Cousins I magnitude error (1)
  99-114  A16   ---     Tel/Inst  Telescope/Instrument
--------------------------------------------------------------------------------
Note (1): Colour-corrected U and S-corrected BVRI Johnson-Cousins magnitudes.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8  F8.2  d       Epoch     Epoch (JD-2400000)
  11- 16  F6.2  d       Phase     Days since explosion
  19- 24  F6.3  mag     umag      ? SDSS u magnitude (2)
  27- 32  F6.3  mag   e_umag      ? SDSS u magnitude error (2)
  35- 40  F6.3  mag     gmag      ? SDSS g magnitude (2)
  43- 48  F6.3  mag   e_gmag      ? SDSS g magnitude error (2)
  51- 56  F6.3  mag     rmag      SDSS r magnitude (2)
  59- 64  F6.3  mag   e_rmag      SDSS r magnitude error (2)
  67- 72  F6.3  mag     imag      SDSS i magnitude (2)
  75- 80  F6.3  mag   e_imag      SDSS i magnitude error (2)
  83- 88  F6.3  mag     zmag      ? SDSS z magnitude (2)
  91- 96  F6.3  mag   e_zmag      ? SDSS z magnitude error (2)
  99-114  A16   ---     Tel/Inst  Telescope/Instrument
--------------------------------------------------------------------------------
Note (2): Colour-corrected u and S-corrected ugriz (unprimed) SDSS magnitudes.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table6.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  8  F8.2  d       Epoch     Epoch (JD-2400000)
  11- 16  F6.2  d       Phase     Days since explosion
  19- 24  F6.3  mag     Jmag      ? 2MASS J magnitude (3)
  27- 32  F6.3  mag   e_Jmag      ? 2MASS J magnitude error (3)
  35- 40  F6.3  mag     Hmag      2MASS H magnitude (3)
  43- 48  F6.3  mag   e_Hmag      2MASS H magnitude error (3)
  51- 56  F6.3  mag     Kmag      ? 2MASS K magnitude (3)
  59- 64  F6.3  mag   e_Kmag      ? 2MASS K magnitude error (3)
  67- 82  A16   ---     Tel/Inst  Telescope/Instrument
--------------------------------------------------------------------------------
Note (3): S-corrected JHK 2MASS magnitudes.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table10.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  8  F8.2  d         Epoch     Epoch (JD-2400000)
  11- 16  F6.2  d         Phase     Days since explosion
  19- 24  F6.3 10+34W     Lum       Pseudo-bolometric UV-MIR (1900-50000{AA})
                                     luminosity (10^41^erg/s)
  27- 32  F6.3 10+34W   e_Lum       Luminosity error
  35- 40  F6.3 10+34W    eLum       Systematic error - Low (4)
  43- 48  F6.3 10+34W    ELum       Systemetic error - High (4)
--------------------------------------------------------------------------------
Note (4): Systematic error arising from the distance and extinction.
--------------------------------------------------------------------------------

Acknowledgements:
    Mattias Ergon, mattias.ergon(at)astro.su.se

================================================================================
(End)  Mattias Ergon [Stockholm University], Patricia Vannier [CDS]  30-Dec-2013
