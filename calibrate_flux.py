from __future__ import print_function
import argparse
import numpy as np


##### Options #######

caltablefile = 'CFHT-K2C9-Microlensing-Calibration_v1.0.pkl'
defaultpos = [1024,2322] #Use this if location on chip not provideded/known



# Process input arguments

parser = argparse.ArgumentParser(description='Calibrate a set of CFHT-K2C9 microlensing survey fluxes.')
parser.add_argument('-e','--event', action='store', dest='event', default=None, help='Event name (e.g., OB160001)')
parser.add_argument('-f','--field', action='store', dest='field', default=None, help='Field name (e.g., CF1)')
parser.add_argument('-c','--chip', action='store', dest='chip', default=None, help='Chip name (e.g., ccd05)')
parser.add_argument('-a','--ra', action='store', dest='ra', default=np.nan, help='RA in degrees',type=float)
parser.add_argument('-d','--dec', action='store', dest='dec', default=np.nan, help='Dec in degrees',type=float)
parser.add_argument('--fluxfile', action='store', dest='fluxlist', default=None, help='File containing multiple lines with ra, dec, field, chip and fluxes.')
parser.add_argument('--save', action='store_true', dest='savefile', help='Save the file as <input>.ps1calf?')
parser.add_argument('gflux', default=None, nargs='?', help='Flux in g filter',type=float)
parser.add_argument('gerror', default=None, nargs='?', help='Error on g flux',type=float)
parser.add_argument('rflux', default=None, nargs='?', help='Flux in r filter',type=float)
parser.add_argument('rerror', default=None, nargs='?', help='Error on r flux',type=float)
parser.add_argument('iflux', default=None, nargs='?', help='Flux in i filter',type=float)
parser.add_argument('ierror', default=None, nargs='?', help='Error on i flux',type=float)
args = parser.parse_args()

fluxdata = np.array([[args.gflux,args.gerror,args.rflux,args.rerror,args.iflux,args.ierror]])
#print(fluxdata.shape)
fluxes = fluxdata[:,0::2]
#print("fluxes",fluxes)
errors = fluxdata[:,1::2]
#print("errors",errors)
fluxespassed=False
field=args.field
chip=args.chip
ra=np.array([args.ra])
dec=np.array([args.dec])
event=args.event
fluxlist=args.fluxlist


import astropy.io.fits as fits
import pickle
import poly2d
import fitsh
import findCFHTevent as finder
import sys
import re

#Figure out if there are incompatable option combos

if not all(x is None for x in fluxes[0]):
    fluxespassed=True

#print(fluxespassed)
#print(event)

if fluxlist==None:

    if fluxespassed==False and event==None:
        parser.print_help()
        sys.exit()
    
    if ((field==None) != (chip==None)) and event==None:
        print("Error: Must provide both field and chip if not specifying event.")
        sys.exit(1)

    if (not np.isfinite(ra[0])) != (not np.isfinite(dec[0])):
        print("Error: Must provide both RA and Dec or neither.")
        sys.exit(1)

if event != None and fluxlist != None:
    print("Error: pass only an event or a fluxlist, not both.")
    sys.exit()

if fluxespassed and fluxlist != None:
    print("Error: pass only individual fluxes or a fluxlist, not both.")
    sys.exit()




eventdata = dict()

#If an event is specified
if event!=None:

    #Event has been specified - we'll use it to determine field, chip, ra, dec, etc.
    eventdata = finder.find_event(event,field)
    checkevent = finder.find_event_errors(eventdata)
    if checkevent != None:
        print(checkevent)
        sys.exit()

    #OK, we have everything that's needed
    event = eventdata['event']
    field = eventdata['field']
    chip = eventdata['chip']
    ra[0] = eventdata['ra']
    dec[0] = eventdata['dec']

    #Look for a file with fluxes if not passed
    if not fluxespassed:
        srcfluxfile = '%s_%s_i/%s/%s.source_flux' % (field,chip,event,event)
        fluxdata = np.loadtxt(srcfluxfile,ndmin=2)

if fluxlist != None:
    fluxdata = np.loadtxt(fluxlist)

#print(fluxdata.shape)
    
#Load the calibration data
caltable = pickle.load(open(caltablefile,'rb'))

#Now we can finally work on the calibrations
instmag=dict()
insterr=dict()
intermag=dict()
calphot=dict()
calerr=dict()
pxy=dict()

fdo=0
if fluxdata.shape[1]==8:
    fdo=2
    ra = fluxdata[:,0]
    dec = fluxdata[:,1]

#print(ra)
#print(dec)

for i,f in enumerate(['g','r','i']):

    fcf = '%s_%s_%s' % (field,chip,f)

    #Pull the data from the pickle
    dophot_zpt,calcolor,g0,g1,flat,astrometry_data,astrometry_dx,astrometry_dy =caltable[fcf]

    #Convert fluxes to instrumental magnitudes
    instmag[f] = dophot_zpt - 2.5*np.log10(fluxdata[:,fdo+2*i])
    insterr[f] = 2.5/np.log(10) * fluxdata[:,fdo+2*i+1]/fluxdata[:,fdo+2*i]

    #print(f,"flux",fluxdata[:,fdo+2*i],"error",fluxdata[:,fdo+2*i+1])

    #Compute the pixel coordinates
    x = np.array([defaultpos[0]])
    y = np.array([defaultpos[1]])

    if ra[0]!=None:
        #print(ra)
        
        pxy[f] = fitsh.ad2pix(ra,dec,astrometry_data,astrometry_dx,astrometry_dy)
        #print(pxy[f])
        x = (pxy[f])[:,0]/poly2d.xscale - 0.5
        y = (pxy[f])[:,1]/poly2d.yscale - 0.5

    #Apply the flatfield correction
    intermag[f] = instmag[f] - poly2d.polyval2d(x,y,flat)    


#Need all the colors to be corrected before applying the color correction
for i,f in enumerate(['g','r','i']):

    fcf = '%s_%s_%s' % (field,chip,f)

    #Pull the data from the pickle
    dophot_zpt,calcolor,g0,g1,flat,astrometry_data,astrometry_dx,astrometry_dy =caltable[fcf]

    c1=calcolor[0]
    c2=calcolor[2]
    color = intermag[c1]-intermag[c2]
    colorerr2 = insterr[c1]**2 + insterr[c2]**2

    calmag = 99 + np.zeros(instmag[f].shape)
    calerror = 99 + np.zeros(instmag[f].shape)
    #Only perform valid computations
    #print(np.isfinite(instmag[c1]))
    mask = np.logical_and(np.isfinite(instmag[c1]),np.isfinite(instmag[c2]))
    mask = np.logical_and(mask,np.isfinite(instmag[f]))

    calmag[mask] = instmag[f][mask] + g0 + g1*color[mask]
    calerror[mask] = np.sqrt(insterr[f][mask]**2 + g1**2*colorerr2[mask])

    calphot[f] = np.copy(calmag)
    calerr[f] = np.copy(calerror)



#Print out the data
fileopen=0
if args.savefile and event!=None and not fluxespassed:
    srcmagfile = '%s_%s_i/%s/%s.source_flux.ps1calf' % (field,chip,event,event)
    fh = open(srcmagfile,'w')
    sys.stdout = fh
    fileopen=1
    

print("#ra dec g g_err r r_err i i_err xi yi i_inst i_inst_err i_type xg yg g_inst g_inst_err g_type xr yr r_inst r_inst_err r_type")
for i,line in enumerate(fluxdata):
    
    print(ra[i],dec[i], \
          (calphot['g'])[i],(calerr['g'])[i], \
          (calphot['r'])[i],(calerr['r'])[i], \
          (calphot['i'])[i],(calerr['i'])[i], \
          (pxy['i'])[i,0],(pxy['i'])[i,1],(instmag['i'])[i],(insterr['i'])[i],-1, \
          (pxy['g'])[i,0],(pxy['g'])[i,1],(instmag['g'])[i],(insterr['g'])[i],-1, \
          (pxy['r'])[i,0],(pxy['r'])[i,1],(instmag['r'])[i],(insterr['r'])[i],-1)


if fileopen==1:
    fh.close()
    sys.stdout = sys.__stdout__
    
