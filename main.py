#!/usr/bin/env python
"""Finding correlation functions using TreeCorr."""
import treecorr
def handle_exception(exc_type, exc_value, exc_traceback):
    """ Handle Exceptions in log."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error("Uncaught exception",
                  exc_info=(exc_type, exc_value, exc_traceback))


if __name__ == "__main__":
    import argparse
    import sys
    import logging
    maskno = 4
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="Set Log Level", type=str,
                        default='INFO')
    args = parser.parse_args()
    loglevel = args.log.upper()
    logging.basicConfig(format='(%(asctime)s %(filename)s %(levelname)s) '
                        + '%(funcName)s %(lineno)d >> %(message)s',
                         filename=f"thelog{maskno}.log",
                        filemode='w',
                        level=getattr(logging, loglevel, None))
    logging.captureWarnings(True)
    sys.excepthook = handle_exception

    rcsconfig = {'ra_col': 'RA',
                 'dec_col': 'DEC',
                 'g1_col': 'e1',
                 'g2_col': 'e2',
                 'w_col': 'w',
                 'ra_units': 'deg',
                 'dec_units': 'deg'}

    szconfig = {'ra_col': 'RA',
                'dec_col': 'DEC',
                'k_col': 'y',
                'ra_units': 'deg',
                'dec_units': 'deg'}

    corrconfig = {'min_sep': 0.25,
                  'max_sep': 1.5,
                  'nbins':10}

    m1config = {'ra_col': 'RA',
                'dec_col': 'DEC',
                'k_col': 'm1',
                'w_col': 'w',
                'ra_units': 'deg',
                'dec_units': 'deg'}


    nside = 2048
    logging.info(f"Calculating for NSIDE = {nside}")
    rcscat = treecorr.Catalog("rcslens.fits", rcsconfig)
    szcat = treecorr.Catalog(f"szmaps_masked{maskno}_{nside}.fits", szconfig)
    m1cat = treecorr.Catalog("rcslens.fits", m1config)

    kg = treecorr.KGCorrelation(corrconfig,logger=logging.getLogger())
    kg.process(szcat, rcscat)   # Calculate the cross-correlation
    kg.write(f"crosscorr{maskno}_{nside}.result")

    nk = treecorr.NKCorrelation(corrconfig,logger=logging.getLogger())
    nk.process(szcat, m1cat)
    nk.write(f"calib{maskno}_{nside}.result")

    ny = treecorr.NKCorrelation(corrconfig,logger=logging.getLogger())
    ny.process(rcscat, szcat)
    ny.write(f"ycorr{maskno}_{nside}.result")
    
    ng = treecorr.NGCorrelation(corrconfig,logger=logging.getLogger())
    ng.process(rcscat, rcscat)
    ng.write(f"shearcorr{maskno}_{nside}.result")

    gg = treecorr.GGCorrelation(corrconfig,logger=logging.getLogger())
    gg.process(rcscat, rcscat)
    gg.write(f"shear_auto_corr{maskno}_{nside}.result")

    kk = treecorr.KKCorrelation(corrconfig,logger=logging.getLogger())
    kk.process(szcat, szcat)
    kk.write(f"sz_auto_corr{maskno}_{nside}.result")

    logging.info("DONE")
