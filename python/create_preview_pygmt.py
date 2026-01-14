#!/usr/bin/env python3
import argparse
import lics_vis as lv
import lics_processing as lp
import os

def main():
    parser=argparse.ArgumentParser(
        description='Creates a preview png for input geotiff'
    )
    parser.add_argument("--grid", required=True,
                        help="input grid to plot (must be a geotiff file, in WGS-84)")
    parser.add_argument("--title", default='preview',
                        help="title")
    parser.add_argument("--label", default='unit',
                        help="label below the colour scale")
    parser.add_argument("--lims", nargs=2, type=float,
                        help="min max (set as separate numbers, e.g. --lims -10 10")
    parser.add_argument("--cmap", default='roma',
                        help="colour scale map (def: roma, try e.g. vik, viridis")
    parser.add_argument("--photobg", action='store_true',
                        help="adds photo background")
    parser.add_argument("--dpi", default=150, type=int,
                        help="DPI for output png")
    parser.add_argument("--projection", default='M10c',
                        help="standard GMT projection string")
    #
    args = parser.parse_args()
    # checks:
    if not args.grid.endswith('tif'):
        parser.error('you must provide a GeoTIFF grid file')
    #
    if not args.lims:
        lims = None
    else:
        lims = args.lims
    ##
    outfile = os.path.splitext(args.grid)[0]+'.png'
    #
    #
    if args.title == 'preview':
        args.title = os.path.splitext(os.path.basename(args.grid))[0]
    #
    ########## generate it
    intif = lp.load_tif2xr(args.grid)
    #
    gg=lv.pygmt_plot(intif,
        title=args.title,
        label=args.label,
        lims=lims,
        cmap=args.cmap,
        projection=args.projection,
        photobg=args.photobg,
        )
    gg.savefig(outfile, dpi=dpi)
    #
    print(outfile)


if __name__ == '__main__':
    main()
