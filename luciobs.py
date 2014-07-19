#! /usr/bin/env python
import pyfits
import numpy as np
from glob import glob
import time


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


def skysub(images, indx1, indx2):
    """Subtract one dither from the previous dither"""
    im1,hdr1 = pyfits.getdata(images[indx1], header=True)
    im2,hdr2 = pyfits.getdata(images[indx2], header=True)
    new = im2 - im1
    
    # get output file name
    base1 = images[indx1].split('.')[2]
    base2 = images[indx2].split('.')[2]
    output = 'im' + base1 + '-' + base2 + '.fits'

    # subtract im2 from im1. Stars from im2 will appear as negative objs
    # output has im1's header
    pyfits.writeto(output, new, header=hdr1, clobber=True)


def main():
    images = glob('luci*.fits')
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

    # 'sky subtract' one image from another
#    indx1 = range(len(images))
#    indx2 = np.roll(indx1, -1)
#
#    # subtract the last image of the night with the one immediately 
#    # preceding it. doesn't make sense to subtract first image from
#    # last image because sky could change over such a long time
#    indx2[-1] = indx2[-3]

    # dither subtract in pairs
    # create array of all even index images
    indx = range(len(images))
    inde = indx[::2]
    # create array of all odd index images
    indo = indx[1::2]
    # if there is an odd number of images, append the last element 
    # of indo to indo so the last image is dither subtracted as well
    indo.append(indo[-1])

    indxpair = zip(inde,indo)
    for id1, id2 in indxpair:
        print id1,id2
        skysub(images, id1, id2)


if __name__ == '__main__':
    main()

