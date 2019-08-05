#!/usr/bin/env python
"""Finding correlation functions using TreeCorr."""
import treecorr
import healpy as hp
import numpy as np
import astropy.io.fits as pf
from read_data import declratoindex
from tqdm import tqdm
source = "/user1/pranav/msc_codes"
NSIDE = 2048


def declratoindex(decl, ra, nside=NSIDE):
    """
    Return the corresponding index in the Healpy array.

    Parameters
    ----------
    decl: Float
            Declination
    RA: Float
            Right Ascention
    Returns
    -------
    index: Int
            The index in the Healpy Array
    Raises
    ------
    None
    See Also
    --------
    None

    """
    assert decl > -90 and decl < 90
    assert ra > 0 and ra < 360
    return hp.pixelfunc.ang2pix(nside, np.radians(90. - decl), np.radians(ra))

def indextodeclra(index, nside=NSIDE):
    """
    Convert index to angles.

    Parameters
    ----------
    index: Int
        Healpix pixel index

    Returns
    -------
    Decl, RA: Float
        Declination and right ascention

    Raises
    ------
    None

    See Also
    --------
    DeclRaToIndex()

    Notes
    -----
    None

    """
    decl, ra = np.degrees(hp.pix2ang(nside, index))
    decl = 90. - decl
    assert np.all(decl >= -90) and np.all(decl <= 90)
    assert np.all(ra >= 0) and np.all(ra <= 360)
    return decl, ra


def make_calib_rcs_fits():
    """ Make RCS Fits file for calibration Purposes"""
    logging.info("Reading RCS Data")
    kk = np.loadtxt(source+"/kids_data/rcslens.csv", delimiter=",",
                    skiprows=1)
    ra = kk[:, 0]
    dec = kk[:, 1]
    e1 = kk[:, 2]
    e2 = kk[:, 3]
    w = kk[:, 4]
    m = kk[:, 5]
    logging.info("Finished Reading Data")
    col1 = pf.Column(name="RA", format='D', array=ra)
    col2 = pf.Column(name="DEC", format='D', array=dec)
    col3 = pf.Column(name="e1", format='D', array=e1)
    col4 = pf.Column(name="e2", format='D', array=e2)
    col5 = pf.Column(name="w", format='D', array=w)
    col6 = pf.Column(name="m", format='D', array=m)
    col7 = pf.Column(name="m1", format='D', array=m+1)
    col8 = pf.Column(name="nege1", format='D', array=-e1)
    col9 = pf.Column(name="nege2", format='D', array=-e1)
    lenspairs = np.stack((e1,e2),axis=-1)
    np.random.shuffle(lenspairs)
    col10 = pf.Column(name="shufe1", format='D', array=lenspairs[:,0])
    col11 = pf.Column(name="shufe2", format='D', array=lenspairs[:,1])
    cols = pf.ColDefs([col1, col2, col3, col4, col5, col6,
                       col7, col8, col9, col10, col11])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto("rcslens_forcalib.fits")
    pf.info("rcslens_forcalib.fits")
    logging.info(pf.open("rcslens_forcalib.fits")[1].header)


def make_rcs_fits():
    """Convert the RCS CSV to FITS File."""
    logging.info("Making RCS Fits file from RCS Data")
    kk = np.loadtxt(source+"/kids_data/rcslens.csv", delimiter=",",
                    skiprows=1)
    ra = kk[:, 0]
    dec = kk[:, 1]
    e1 = kk[:, 2]
    e2 = kk[:, 3]
    w = kk[:, 4]
    m = kk[:, 5]
    logging.info("Finished Reading Data")

    """
    szmap = hp.read_map(source+'/kids_data/COM_CompMap_Compton-SZMap-milca-ymaps_2048_R2.00.fits')
    szlist = []
    for i in tqdm(range(len(ra))):
        index = int(declratoindex(dec[i], ra[i], NSIDE))
        szlist.append(szmap[index])
    szlist = np.array(szlist)

    shearcalib = np.ones(m.shape) + m
    e1_calib = e1/shearcalib
    e2_calib = e2/shearcalib
    w_calib = w*shearcalib
    """
    col1 = pf.Column(name="RA", format='D', array=ra)
    col2 = pf.Column(name="DEC", format='D', array=dec)
    col3 = pf.Column(name="e1", format='D', array=e1)
    col4 = pf.Column(name="e2", format='D', array=e2)
    col5 = pf.Column(name="w", format='D', array=w)
    col6 = pf.Column(name="m", format='D', array=m)
    col7 = pf.Column(name="m1", format='D', array=m+1)
    col8 = pf.Column(name="nege1", format='D', array=-e1)
    cols = pf.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto("rcslens.fits")
    pf.info("rcslens.fits")
    logging.info(pf.open("rcslens.fits")[1].header)


def make_sz_map_fits(nside=NSIDE):
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    logging.info(f"For NSIDE = {nside}")
    szmap = hp.read_map(source+'/kids_data/nilc_ymaps.fits') #Just read maps reads the first column. 
                                                             # Verified that it is the same as hdu[1].data['FULL']

    maskfile = pf.open(source+'/kids_data/planck_mask.fits')
    maskno = 4
    mask = maskfile[1].data[f"M{maskno}"]
    szmap = hp.ud_grade(szmap,nside)
    bmap = []
    lmap = []
    ymap = []
    for i in tqdm(range(len(szmap))):
        if mask[i]:
            b, l = indextodeclra(i,nside=nside)
            lmap.append(l)
            bmap.append(b)
            ymap.append(szmap[i])

    lmap = np.array(lmap)
    bmap = np.array(bmap)
    gc_sz = SkyCoord(l=lmap, b=bmap, frame='galactic', unit='deg')
    gcout_sz=gc_sz.transform_to('icrs',merge_attributes=False)
    declmap, ramap = gcout_sz.dec.degree, gcout_sz.ra.degree

    ymap = np.array(ymap)
    print(ymap)
    print("Min",np.min(ymap))
    print("Max",np.max(ymap))
    print("Mean",np.mean(ymap))
    print("Std",np.std(ymap))
    exit()
    assert len(ymap) == len(ramap) == len(declmap)
    logging.info(f"Total Shape : {ymap.shape}")
    col1 = pf.Column(name="RA", format='D', array=ramap)
    col2 = pf.Column(name="DEC", format='D', array=declmap)
    col3 = pf.Column(name="y", format='D', array=ymap)
    cols = pf.ColDefs([col1, col2, col3])
    tbhdu = pf.BinTableHDU.from_columns(cols)
    tbhdu.writeto(f"szmaps_masked{maskno}_eq_{nside}.fits")
    pf.info(f"szmaps_masked{maskno}_eq_{nside}.fits")
    logging.info(pf.open(f"szmaps_masked{maskno}_eq_{nside}.fits")[1].header)


def gal2eqcoords(b, l):
    logging.info("Converting Coordinate Systems")
    try:
        assert (l >= -90).all()
    except AssertionError:
        logging.info(f"{min(l)} > {-90}")
        exit()
    try:
        assert (l <= 90).all()
    except AssertionError:
        logging.info(f"{max(l)} < {90}")
        exit()
    try:
        assert (b >= 0).all()
    except AssertionError:
        logging.info(f"{min(b)} > 0")
        exit()
    try:
        assert (b <= 360).all()
    except AssertionError:
        logging.info(f"{max(b)} < {360}")
        exit()

    b = np.radians(b)
    l = np.radians(l)
    declG = np.radians(np.full(b.shape, 27.12825))
    alphaG = np.radians(np.full(b.shape, 192.85948))
    l_NCP = np.radians(np.full(b.shape, 122.93192))
    sin_decl = np.sin(b)*np.sin(declG) + np.cos(b)*np.cos(declG)*np.cos(l_NCP - l)
    decl = np.arcsin(sin_decl)

    cos_ra = (np.sin(b)*np.cos(declG) - np.cos(b)*np.sin(declG)*np.cos(l_NCP - l))/np.cos(decl)

    sin_ra = np.cos(b)*np.sin(l_NCP - l)/np.cos(decl)
    
    ra = np.arctan2(sin_ra, cos_ra) + alphaG 
    ra = ra - min(ra)
    # decl = decl - min(decl) - (np.pi/2)
    try:
        assert (np.pi/2 >= decl).all()
    except AssertionError:
        logging.info(f"{max(decl)} > {np.pi/2}")
        logging.info(ra[decl>np.pi/2])
        exit()
    try:
        assert (-np.pi/2 <= decl).all()
    except AssertionError:
        logging.info(f"{min(decl)}<{-np.pi/2}")
        logging.info(ra[decl< -np.pi/2])
        exit()
    try:
        assert (2*np.pi >= ra ).all()
    except AssertionError:
        logging.info(f"{max(ra)} > {2*np.pi}")
        logging.info(ra[ra > 2*np.pi])
        exit()
    try:
        assert (0 <= ra).all()
    except AssertionError:
        logging.info(f"{min(ra)} < 0")
        logging.info(ra[ra < 0])
        exit()
    logging.info(f"All Checks Done")
    return np.degrees(decl), np.degrees(ra)
    
    

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="Set Log Level", type=str,
                        default='INFO')
    args = parser.parse_args()
    loglevel = args.log.upper()
    logging.basicConfig(format='(%(asctime)s %(filename)s %(levelname)s) '
                        + '%(funcName)s %(lineno)d >> %(message)s',
                        # filename=f"thelog_makefits.log",
                        # filemode='w',
                        level=getattr(logging, loglevel, None))
    logging.captureWarnings(True)
    sys.excepthook = handle_exception
    #for nside in [2048, 1024, 512, 256]:
        #make_sz_map_fits(nside)
    # make_rcs_fits()
    make_calib_rcs_fits()
    # make_sz_map_fits()
