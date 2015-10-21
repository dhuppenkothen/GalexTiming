import argparse

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from gPhoton.gFind import gFind
from gPhoton.gAperture import gAperture

def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text

def make_int(text):
    return int(text.strip('" '))


def read_catalogue(filename, band="NUV", dt=10.0, target_dir="./", start_ind=0, end_ind=None):
    """
    Read Galex Data for some of the catalogues I've got.
    Note: This will only work with the LMXB and HMXB catalogue files in the 
    data directory!
    
    Everything else is probably not going to be in the right format. However, 
    the code below should be adaptable to whatever format your catalogue is in!
    
    Parameters
    ----------
    filename: string
        The file name (including path) to your catalogue file
    
    band: {"NUV | "FUV"}, optional, default "NUV"
        The band to scrape the data from; either NUV (near UV) or FUV (far UV)
        TODO: check whether there are more than two bands!
    
    dt: float, optional, default 10
        The desired time resolution of the output light curve
    
    target_dir: string, optional, default "./"
        The target directory for the output data, default is the current working directory
    
    start_ind: int, optional, default 0
        If you want to load fewer rows from the catalogue file, this is the 
        index of the row you'd like to start with
    end_ind: int, optional, default None
        If you want to load fewer rows from the catalogue file, this is the 
        index of the row you'd like to end with.
        If None, then this will be set to the last row in the catalogue.
        
    """
    
    ## read in the catalogue
    cat = pd.read_csv(filename, sep="|", skipinitialspace=True,
                        names=["name", "ra", "dec", "vmag", "bv_color", "porb", "flux_limit",
                              "flux", "flux_max", "xray_type", "pulse_period", "alt_name_1",
                              "alt_name_2"],
                        converters = {'name' : strip,
                                      'ra' : strip,
                                      'dec' : strip,
                                      'vmag' : strip,
                                      'bv_color': strip,
                                      "porb": strip,
                                      "flux_limit": strip,
                                      "flux": strip,
                                      "flux_max": strip,
                                      "xray_type": strip,
                                      "pulse_period": strip,
                                      "alt_name_1": strip,
                                      "alt_name_2": strip},
                       usecols=range(1,14,1), index_col=False,
                       skiprows=1)
    
    ## extract RA and Dec:
    ra_all = np.array(cat["ra"])
    dec_all = np.array(cat["dec"])
    
    ## store in coords object:
    coords_all = [SkyCoord("%s %s"%(ra, dec), unit=(u.hourangle, u.deg)) \
              for ra,dec in zip(ra_all, dec_all)]
    

    ## if the last index is None, then the last item of the 
    ## coordinate list is the final index to search through
    if end_ind is None:
        end_ind = len(coods_all)
        
        
    for i,c in enumerate(coords_all[start_ind:end_ind]):
        ## remove white spaces from object identifier
        obj_id = "".join(cat.loc[i,"name"].split())
        search_galex(c, band, dt, obj_id, target_dir)
    return

def search_galex(coords, band, dt, obj_id, target_dir="./"):
    """
    Search the GALEX gPhoton database for observations at a specific set 
    of coordinates, and if data exists, download it and save it to disc.
    
    Parameters
    ----------
    coords: astropy.coordinates.SkyCoord object
        an object with the sky coordinates where to look for GALEX data
    
    band: {"NUV" | "FUV"}
        Look for data in the near-UV or far-UV?
    
    dt: float
        The desired time resolution of the output light curve
    
    obj_id: string
        The object identifier to be used to save the data
        
    target_dir: string, optional, default "./"
        The directory to store the data in. By default the current working directory
    
    """
    ## use gFind to see if there's data in GALEX
    res = gFind(band=band,skypos=[coords.ra.degree,coords.dec.degree])
    ## if not, continue
    if res[band]["expt"] == 0:
        return
    else:
        tstart = res[band]["t0"]
        tend = res[band]["t1"]
        #print(tstart[0])
        #print(tend[0])
        for i,(ts, te) in enumerate(zip(tstart, tend)):
            ## otherwise find the object's name and fetch the data
            d = gAperture(band=band, skypos=[coords.ra.degree,coords.dec.degree],radius=0.03,
                          annulus=[0.03,0.04], trange=[ts, te], stepsz=dt,
                          csvfile="%s%s_%s_lc%i.dat"%(target_dir, obj_id, band, i))
            ## move photons into a DataFrame
            df = pd.DataFrame(d["photons"])
            ## store in HDF5 file:
            df.to_hdf('%s%s_%s_%i_photons.h5'%(target_dir,obj_id, band, i),
                      'table',append=False)

    return

def parse_args():
    parser = argparse.ArgumentParser(description="Read data from a catalogue")

    parser.add_argument("-f", "--filename", action="store", dest="filename", type=str,
                        required=True, help="Catalogue file name and directory")

    parser.add_argument("-b", "--band", action="store", dest="band", type=str,
                        required=True, help="Either 'NUV' (near-UV) or 'FUV' (far-UV) band")
 
    parser.add_argument("-t", "--dt", action="store", dest="dt", type=float,
                        required=False, default=10.0, help="The time resolution of the light "+ 
                        "curve (default 10 seconds)")
 
    parser.add_argument("-d", "--targetdir", action="store", dest="target_dir", type=str,
                       required=False, default="./", help="The target directory for the data (default: PWD)")
 
    parser.add_argument('-s', '--startind', action="store", dest="start_ind", type=int,
                        required=False, default=0, help="The index of the row in the catalogue to start with")

    parser.add_argument("-e", "--endind", action="store", dest="end_ind", required=False, default="None",
                        help="The index of the row in the catalogue to end with (default: last row")

    clargs = parser.parse_args()

    if clargs.end_ind == "None":
        clargs.end_ind = None
    else: 
        clargs.end_ind = int(clargs.end_ind)
    return clargs



if __name__ == "__main__":
    clargs = parse_args()

    print(clargs)

    read_catalogue(clargs.filename, band=clargs.band, dt=clargs.dt, target_dir=clargs.target_dir, 
                   start_ind=clargs.start_ind, end_ind=clargs.end_ind)

