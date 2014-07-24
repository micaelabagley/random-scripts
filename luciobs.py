#! /usr/bin/env python
import pyfits
import numpy as np
from glob import glob
import time
from align import align
import pyregion
from astropy.table import Table
from aperphot import apphot
from astropy.modeling import models,fitting
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def flatcombine(ims):
    """median combine a set of science images of different dithers 
    into a master flat"""
    output = 'Flat.fits'
    nflats = len(ims)
    if nflats < 3:
        print 'There are not enough flats to median combine (fewer than 3)'
        exit()

    # get basic header info
    hdr = pyfits.getheader(ims[0])
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    median_array = np.zeros((nflats,ny,nx), dtype=float)

    for i,image in enumerate(ims):
        im,ihd = pyfits.getdata(image, header=True)
        med = np.median(im)
        if med == 0:
            print 'Median value for ' + im + 'is 0!'
            exit()

        # normalize by median value
        median_array[i,:,:] = im / med

    flat = np.median(median_array, axis=0)

    # add keywords to header
    hdr['NCOMBINE'] = (nflats, 'Number of flats combined to form master flat')
    for i in range(nflats):
        hdr['IMLIST%i'%(i+1)] = (ims[i], 'Image %i used to create flat'%(i+1))

    pyfits.writeto(output, flat, header=hdr, clobber=True)


def flatfield(imlist, flat):
    """Flat-field all images"""
    # check the flat for zeros, replace with median of flat
    flat[flat == 0] = np.median(flat)
    flat = flat / np.median(flat)

    for image in imlist:
        im,hdr = pyfits.getdata(image, header=True)

        # check to see if image has already been flat-fielded
        test = hdr.get('FLATCOR', '')
        if len(test) > 0:
            print '   ' + image + ' already flat-fielded. Skipping '
        elif len(test) == 0:
            # get the date and time
            now = time.strftime('%c')

            # write a new keyword to flag image as having been flat-fielded
            flatstr = 'flat-fielded: ' + now
            hdr['FLATCOR'] = flatstr
            
            new = im / flat

            output = image.split('.fits')[0] + '.ff.fits'
            pyfits.writeto(output, new, header=hdr, clobber=True)


def expt_norm(imlist):
    for image in imlist:
        im,hdr = pyfits.getdata(image, header=True)
        exptime = hdr['EXPTIME']
        new = im / exptime
        # add a header keyword
        hdr['EXPTNORM'] = (exptime, 'Image normalized by EXPTIME')
        pyfits.writeto(image, new, header=hdr, clobber=True)


def dithersub(image1, image2):
    """Subtract one dither from the previous dither"""
    im1,hdr1 = pyfits.getdata(image1, header=True)
    im2,hdr2 = pyfits.getdata(image2, header=True)
    new1 = im2 - im1
    new2 = im1 - im2

    # trim edges 
    new1 = new1[:1950, :1956]
    new2 = new2[:1950, :1956]
    
    # get output file name
    #base1 = images[indx1].split('.')[2]
    #base2 = images[indx2].split('.')[2]
    #output = 'im' + base1 + '-' + base2 + '.fits'

    # subtract im2 from im1. Stars from im2 will appear as negative objs
    # output has im1's header
    pyfits.writeto('L1544_J_dither1.fits', new1, header=hdr1, clobber=True)
    pyfits.writeto('L1544_J_dither2.fits', new2, header=hdr2, clobber=True)

def reducedata():
    """
    # all images except 066, which has a satelite trail
    images = [i for i in glob('luci*.fits') if i not in \
        ['luci.20130320.0066.fits']]
    images.sort()

    # attempt to make flat from images
    flatlist = [images[0], images[1], images[2], images[4], images[8], \
                images[10], images[11]]
    flatcombine(flatlist)

    # flat field all images
    flat, fhd = pyfits.getdata('Flat.fits', header=True)
    flatfield(images, flat)

    images = glob('luci*.ff.fits')
    images.sort()

    # normalize by exposure time to get counts/s
    expt_norm(images)

    # align all images to image1
    imlist = images[1:]
    align(images[0], imlist)
    # assign original headers to aligned images
    aligned = glob('alipy_out/luci*.fits')
    for al in aligned:
        # get header from original image
        alfile = al.split('alipy_out/')[1]
        orig = ''.join([alfile.split('_affineremap')[0], 
            alfile.split('_affineremap')[1]])
        hdr = pyfits.getheader(orig)
        pyfits.update(al, pyfits.getdata(al), header=hdr)

    raw_input('Change alipy_out to alipy_dither1')

    # align all images to image2 for dither subtraction
    imlist = [im for im in images if im not in [images[1]]]
    align(images[1], imlist)

    # assign original headers to aligned images
    aligned = glob('alipy_out/luci*.fits')
    for al in aligned:
        # get header from original image
        alfile = al.split('alipy_out/')[1]
        orig = ''.join([alfile.split('_affineremap')[0], 
            alfile.split('_affineremap')[1]])
        hdr = pyfits.getheader(orig)
        pyfits.update(al, pyfits.getdata(al), header=hdr)

    raw_input('Change alipy_out to alipy_dither2')

    print 'Combine the images with imcombine.'
    raw_input('Press ENTER to continue...') 
    """
    # dither subtract the combined odds from the combined evens
    dithersub('dither1.fits', 'dither2.fits')


def calc_maglim(im, nx1, nx2):
    """Calculate the limiting magnitude for a given aperture size"""
    # generate random coordinates for apertures
    xx = np.random.rand(100) * nx1
    yy = np.random.rand(100) * nx2

    names = ('xc', 'yc', 'flux', 'eflux')
    dtypes = ('float', 'float', 'float', 'float')
    fill = np.zeros(xx.shape, dtype=float)
    t = Table(data=[xx,yy,fill,fill], names=names, dtype=dtypes)

    # perform aperture photometry
    # 30 pixels in radius and set sky to 0
    table = apphot(im, t, 15, [20,30], setsky=0.)

    # bin the photometry
    n,bins = np.histogram(table['flux'], bins=20)
    # get bin centers 
    bincenters = 0.5*(bins[1:]+bins[:-1])
    width = 1.01 * (bins[1] - bins[0])
    # fit a gaussian
    gaussfunc = lambda x,a,b,c: a * np.exp(-(x-b)**2/(2*c**2))
    p0 = [np.max(n),0.,1.]
    popt,pcov = curve_fit(gaussfunc, bincenters, n, p0=p0)
    
    plt.plot(bincenters, n)
    plt.plot(bincenters, gaussfunc(bincenters, *popt))
    plt.savefig('test.pdf')
    print popt


def gauss_centroid(im, x_init, y_init):
    """Find centroid of source by fitting a 2D Gaussian"""
    fitter = fitting.LevMarLSQFitter()

    size = 40
    # get initial guesses for the parameters 
    amplitude = np.max( im[y_init-size:y_init+size,
                           x_init-size:x_init+size] )
    # super rough guess at width of source
    sig_x = 10.
    sig_y = 10.
    g_init = models.Gaussian2D(amplitude, x_init, y_init, sig_x, sig_y)

    # get indices of subsection of image
    x,y = np.mgrid[x_init-size:x_init+size, y_init-size:y_init+size]

    print x_init, y_init
    print x.shape
    print y.shape

    p = fitter(g_init, x, y, im[y_init-size:y_init+size, 
        x_init-size:x_init+size] )

    cenx = p.x_mean.value
    ceny = p.y_mean.value
    print cenx,ceny
    return cenx,ceny


def photometry(dither=1):
    # read in region file
    if dither == 2:
        image,hdr = pyfits.getdata('L1544_J_dither2.fits', header=True)
        regpix = pyregion.open('L1544_pix_dither2.reg')
    else:
        image,hdr = pyfits.getdata('L1544_J_dither1.fits', header=True)
        regpix = pyregion.open('L1544_pix_dither1.reg')
        nx1 = hdr['NAXIS1']
        nx2 = hdr['NAXIS2']
#        calc_maglim(image, nx1, nx2)

    regsky = pyregion.open('L1544_sky.reg')
    
    names = ('star', 'ra', 'dec', 'xc', 'yc', 'flux', 'eflux', 
             'mag', 'emag')
    dtypes = ('int', 'S11', 'S11', 'float', 'float', 'float', 'float',
              'float', 'float')
    t = Table(data=None, names=names, dtype=dtypes)

    # centroid sources
    cenx = np.zeros(len(regpix), dtype=float)
    ceny = np.zeros(len(regpix), dtype=float)
    for i,star in enumerate(regpix):
        id = star.comment
        ra = regsky[i].coord_list[0]
        dec = regsky[i].coord_list[1]

#        x_init = star.coord_list[0]
#        y_init = star.coord_list[1]
        cenx = star.coord_list[0]
        ceny = star.coord_list[1]
#        # centroid by fitting a 2D Gaussian
#        cenx,ceny = gauss_centroid(image, x_init, y_init)
        
        t.add_row([id, ra, dec, cenx, ceny, 0.0, 0.0, 0.0, 0.0])
    
    # perform aperture photometry
    # 30 pixels in radius and sky annuli at 30-50 pixels
    table = apphot(image, t, 15, [20,30], return_mag=True)

    formats = {'xc':'%8.3f', 'xc':'%8.3f', 'flux':'%7.2f', \
               'eflux':'%6.2f', 'mag':'%6.3f', 'emag':'%5.3f'}

    if dither == 2:
        table.write('flux_dither2.txt', format='ascii.tab', formats=formats)
    else:
        table.write('flux_dither1.txt', format='ascii.tab', formats=formats)


def main():
    #reducedata()
    photometry()
    photometry(dither=2)



if __name__ == '__main__':
    main()

