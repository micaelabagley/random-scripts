#! /usr/bin/env python
#######################################################################
## usage: align.py [-h] ref_file unk_files
##
## Align a list of images to a reference image
##
## positional arguments:
##   ref_file    File name of reference image
##     unk_files   List of images to be aligned
##
##     optional arguments:
##       -h, --help  show this help message and exit
#######################################################################
import argparse
import os
import alipy
import pyfits


def align(ref_file, unk_files):
    """Align images using alipy"""
    # generate list of output files with altered headers
    outfiles = []
    headers = []
    for f in unk_files:
        outroot,outext = os.path.splitext(f)     # peel off extension
        outfile = ''.join([outroot,'_d',outext]) # add _d to filename
        outfiles.append(outfile)

        header = pyfits.getheader(f)
        header['FILENAME'] = os.path.basename(outfile)
        header['ALIGNED'] = (os.path.basename(ref_file),'Reference file')
        headers.append(header)

    # Run alignment on file list
    ident = alipy.ident.run(ref_file, unk_files, visu=False, sexrerun=True)

    # For each unk file, align and make new files
    for idx,imid in enumerate(ident):
        ref = imid.ref
        unk = imid.ukn

        #transform = alipy.star.fitstars(unk.starlist,ref.starlist)
        outputshape = alipy.align.shape(ref_file)
        print unk.filepath
        print outfiles[idx]

        # Shift and write new file
        alipy.align.affineremap(unk.filepath, transform=imid.trans,
            shape=outputshape)#, alifilepath=outfiles[idx])

        # Update new header
        #pyfits.update(outfiles[idx], pyfits.getdata(outfiles[idx]),
        #   header=headers[idx])


def main():
    parser = argparse.ArgumentParser(
        description='Align  a list of images to a reference image')
    parser.add_argument('ref_file', type=str, nargs=1,
        help='File name of reference image')
    parser.add_argument('unk_files', type=str, nargs=1,
        help='List of images to be aligned')
    args = parser.parse_args()

    ref_file = args.ref_file[0]
    unk_files = args.unk_files[0]

    print ref_file
    print unk_files
    exit()
    align(ref_file, unk_files)


if __name__ == '__main__':
    main()

