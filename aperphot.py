#! /usr/bin/env python
import numpy as np
import pyfits
from astropy.table import Table
import matplotlib.pyplot as plt


def format_output():
    names = ('star', 'ra', 'dec', 'xc', 'yc', 'flux', 'eflux')
    dtypes = ('int', 'S11', 'S11', 'float', 'float', 'float', 'float')
    t = Table(data=None, names=names, dtype=dtypes)

    

def aperture(im, xc, yc, radius):
    """Identify pixels inside circular aperture"""
    # get x and y vectors of image
    y,x = np.ogrid[:im.shape[0], :im.shape[1]]
    # distance of every pixel from the star's position
    r2 = (x - xc)**2 + (y - yc)**2
    # find all pixels inside a circular aperture
    mask = r2 <= (radius - 0.5)**2

    # the mast will be applied to the image such that all pixels inside
    # the circle have weight 1 and all outside have weight 0.
    # so ecast as floats
    #return mask.astype(float)
    return mask


def sky_annulus(im, xc, yc, skyrad):
    """Identify pixels inside sky annulus"""
    y,x = np.ogrid[:im.shape[0], :im.shape[1]]
    r2 = (x - xc)**2 + (y - yc)**2
    mask = (r2 <= (skyrad[1]-0.5)**2) & (r2 >= (skyrad[0]-0.5)**2)
    #return mask.astype(float), npix
    return mask


def apphot(im, table, radius, skyrad, setsky=None, return_mag=False):
    """Perform aperture photometry"""

    for i in range(len(table)):
        if i % 10 == 0:
            print i
        xc = table['xc'][i]
        yc = table['yc'][i]
        aper = aperture(im, xc, yc, radius)
        sky = sky_annulus(im, xc, yc, skyrad)

        #plt.imshow(aper,origin='lower')
        apermask = np.ma.MaskedArray(im,~aper)
        counts = np.ma.sum(apermask)
        #ap_area = np.pi * radius**2
        ap_area = apermask.count()

        skymask = np.ma.MaskedArray(im,~sky)

        med_sky = np.ma.median(skymask)
        #print counts,med_sky, counts-(med_sky*ap_area)
        if setsky:
            flux = counts - (setsky * ap_area)
        else:
            flux = counts - (med_sky * ap_area)
        

        #counts = np.sum(aper * im)
        #med_sky = np.median(sky * im)
        #flux = counts - (med_sky * ap_area)
        skyvar = np.std(skymask)**2
        npix = skymask.count()
        sigsq = skyvar / npix

        error1 = ap_area * skyvar # scatter in sky values
        error2 = flux # random photon noise
        error3 = sigsq * ap_area**2 # uncertainty in mean sky brightness
        error = np.sqrt(error1 + error2 + error3)
        if return_mag:
            mag = 24.519 - 2.5*np.log10(flux)
            emag = 1.0857 * error / flux

        # update flux and error in table
        table['flux'][i] = flux
        table['eflux'][i] = error
        if return_mag:
            table['mag'][i] = mag
            table['emag'][i] = emag

    return table



def main():
    apphot()


if __name__ == '__main__':
    main()

