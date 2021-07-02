import numpy as np
import string
import pickle
import os
import sys
import options
import math

#______________________________________________________________________________

def allFields(year):

# 2015 has 18 fields (72 ccds)
  if year=='2015':
    fields=['11','12','13','14',
            '21','22','23','24',
            '31','32','33','34',
            '41','42','43','44',
            '5','6']

# 2016 has 32 fields (128 ccds)
  elif year=='2016':
    fields=['11','12','13','14',
            '21','22','23','24',
            '31','32','33','34',
            '41','42','43','44',
            '51','52','53','54',
            '61','62','63','64',
            '71','72','73','74',
            '81','82','83','84']

# 2017 has 56 fields (224 ccds)
  elif year=='2017':
    fields=['n1_1','n1_2','n1_3','n1_4',
            'n2_1','n2_2','n2_3','n2_4',
            'n3_1','n3_2','n3_3','n3_4',
            'n4_1','n4_2','n4_3','n4_4',
            'c1_1','c1_2','c1_3','c1_4',
            'c2_1','c2_2','c2_3','c2_4',
            'c3_1','c3_2','c3_3','c3_4',
            'c4_1','c4_2','c4_3','c4_4',
            's1_1','s1_2','s1_3','s1_4',
            's2_1','s2_2','s2_3','s2_4',
            's3_1','s3_2','s3_3','s3_4',
            's4_1','s4_2','s4_3','s4_4',
            's5_1','s5_2','s5_3','s5_4',
            's6_1','s6_2','s6_3','s6_4']

# 2018 has 56 fields (224 ccds)
#   (but only 160 have K band.  64 north ones do not have it)
  elif year=='2018':
    fields=['n1_1','n1_2','n1_3','n1_4',
            'n2_1','n2_2','n2_3','n2_4',
            'n3_1','n3_2','n3_3','n3_4',
            'n4_1','n4_2','n4_3','n4_4',
            'c1_1','c1_2','c1_3','c1_4',
            'c2_1','c2_2','c2_3','c2_4',
            'c3_1','c3_2','c3_3','c3_4',
            'c4_1','c4_2','c4_3','c4_4',
            's1_1','s1_2','s1_3','s1_4',
            's2_1','s2_2','s2_3','s2_4',
            's3_1','s3_2','s3_3','s3_4',
            's4_1','s4_2','s4_3','s4_4',
            's5_1','s5_2','s5_3','s5_4',
            's6_1','s6_2','s6_3','s6_4']
  else:
    exit('ERROR: no fields specified for this year')

  ccds=['1','2','3','4']

  return fields,ccds

def getFieldfromName(name):

  year = name[2:6]

  if year=='2015' or year=='2016':
    # actually this doesn't work for 2015 fields 5,6 which only have 1 digit
    #field = name[7:9]
    #ccd = name[10:11]
    field = name.split('_')[1]
    ccd = name.split('_')[2]
  else:
    field = name[7:11]
    ccd = name[12:13]

  return year,field,ccd

def wavebands(year,field):

  if year=='2015' or year=='2016':
    band='H'
    altband=''
  else:
    if 'c' in field or 's1' in field or 's2' in field:
      if 'n' in field or 's3' in field or 's4' in field \
         or 's5' in field or 's6' in field:
        print(' field:',field)
        exit ('ERROR: trouble setting the bandpass1')
      band='K'
      altband='H'
    elif 'n' in field or 's3' in field or 's4' in field \
         or 's5' in field or 's6' in field:
      if 'c' in field or 's1' in field or 's2' in field:
        print(' field:',field)
        exit ('ERROR: trouble setting the bandpass2')
      band='H'
      altband='K'
    else:
      print(' year,field:',year,field)
      exit('ERROR: trouble setting the bandpass3')

  return band,altband
      
def convertRADec(ra_hms,dec_hms):
#  print 'ra dec inputs:',ra_hms,dec_hms
  ra_hms=string.strip(ra_hms)
  ra_hr=float(ra_hms[:2])
  ra_min=float(ra_hms[3:5])
  ra_sec=float(ra_hms[6:])
#  print 'ra',ra_hr,'x',ra_min,'x',ra_sec
  ra=ra_hr+ (ra_min + (ra_sec/60.))/60.
  ra=ra/24.*360.

  dec_hms=string.strip(dec_hms)
  dec_deg=float(dec_hms[:3])
  dec_min=float(dec_hms[4:6])
  dec_sec=float(dec_hms[7:])
#  print 'dec',dec_deg,'x',dec_min,'x',dec_sec
  dec=abs(dec_deg)+ (dec_min + (dec_sec/60.))/60.
  dec=dec*np.sign(dec_deg)
  return ra,dec

def getOGLEdata(year,options,targetSelection='All',
                verbose=False):
  '''
  Read in the photometry from the /data directory.
  If it is not already downloaded, run downloadData.py
  The year to be read in is specified.
  Optionally select a batch of targets; otherwise all targets considered.
  '''
  import os
  from math import log10

  if targetSelection==[]:
    targetSelection='All'
  
#  dir='data/ogle/'+year+'/'
  dir=options['dataDir']+'/ogle/'+year+'/'

  # open the list of all the events and their OGLE single-lens fit params
  filename='lenses.par'
  if not os.access(dir+filename, os.R_OK):
    print(dir+filename)
    exit('ERROR: no OGLE data for this year.  Run downloadData.py')
  file=open(dir+filename,'r')

  events={} 
  for line in file:
    if line.strip()=='':
#      print 'skipping a blank line'
      pass
    else:
      fields=line.split()
      if fields[0]=='Event':
        header=fields
      else:
        name=fields[header.index('Event')]
# an optional list can be used to select specific targets within year
        if targetSelection=='All' or name in targetSelection:
          event={'data':{}}
#                 't0':[],'tE':[],'u0':[],
#                 'fb':[],'Ftot':[],'Fs':[]}
          event['t0']=float(fields[header.index('Tmax(HJD)')])
          event['tE']=float(fields[header.index('tau')])
          event['u0']=float(fields[header.index('umin')])
          if event['u0']==0:
#            print 'NOTE: u0 is zero.  using 1/A instead'
            event['u0']=1./float(fields[header.index('Amax')])
          event['fb']=float(fields[header.index('fbl')])
          event['Ftot']=float(fields[header.index('I_bl')])
          event['Fs']=float(fields[header.index('I0')])
#          print 'name:',event['name']

          if 0:
# verify consistency with redundant info in the table
            Amax=float(fields[header.index('Amax')])
            Dmag=float(fields[header.index('Dmag')])
            u0=event['u0']
#            print 'Amax=1/u0?',name,Amax-0.,1./u0,(Amax-0.)*u0  # not good
#            print 'Dmag=Amax?',name,Amax,Dmag,Amax/10.**(Dmag/2.5) # good
            I0=event['Fs']
            I0calc=event['Ftot'] - 2.5*log10(event['fb'])
#            print 'I0=fbl*Ibl? ',I0,I0calc,I0/I0calc            # good

# now load in the photometry        
          subdir=dir+name+'/'
          if not os.access(subdir+'phot.dat', os.R_OK):
            exit('ERROR: no OGLE data for this target.  Run downloadData.py')
          photomFile=open(subdir+'phot.dat','r')

          if verbose:
            print('loading in the OGLE info:',name)

          photom={'time':[],
                  'mag':[],'magerr':[],'flux':[],'fluxerr':[]}
# use this factor to convert a mag err to flux err
          errorConvertFactor=0.92
# set the I band zeropoint, to convert magnitudes to Jy
          zeropointI=1.
          for line in photomFile:
            fields=line.split()
            if float(fields[2])>0.5:
              if verbose:
                print( ' skipping a bad OGLE flux', \
                      fields[1],fields[2])
            else:
              photom['time'].append(float(fields[0]) -
                                    options['JDadjust'])
              photom['mag'].append(float(fields[1]))
              photom['magerr'].append(float(fields[2]))
          photom['time']=np.array(photom['time'])
          photom['mag']=np.array(photom['mag'])
          photom['magerr']=np.array(photom['magerr'])
            
          photom['flux']=zeropointI* \
                        10.**(-photom['mag']/2.5)
          photom['fluxerr']=photom['magerr']* \
                           photom['flux']*errorConvertFactor

#          photom['Nphotom']=len(photom['time'])
#          print 'Nphotom',ogle['Nphotom'],len(ogle['fluxerr'])
          event['data']['ogle']=photom

          events[name]=event
                  
  print('# of events for',year+':',len(events.keys()))

  #print events
  
# if specific targets are requested, verify that they were found  
  #print 'targetSelection',targetSelection
  if targetSelection!='All':
    for targetName in targetSelection:
      if not targetName in events.keys():
        if targetName!='NOGLE':
          print('TROUBLE: selected target not found!',targetName)
          exit('blank data will crash later!')

  return events
#______________________________________________________________________________

def getOGLEheader(year,options,targetSelection='All'):
  '''
  Read in just the header info for OGLE alerts
  The year to be read in is specified.
  Optionally select a batch of targets; otherwise all targets considered.
  '''
  import os

  if targetSelection==[]:
    targetSelection='All'
  
#  dir='data/ogle/'+year+'/'
  dir=options['dataDir']+'/ogle/'+year+'/'
# open the list of all the events and their OGLE single-lens fit params
  filename='lenses.par'
  if not os.access(dir+filename, os.R_OK):
    exit('ERROR: no OGLE data for this year.  Run downloadData.py')
  file=open(dir+filename,'r')

  events={} 
  for line in file:
    if line.strip()=='':
#      print 'skipping a blank line'
      pass
    else:
      fields=line.split()
      if fields[0]=='Event':
        header=fields
      else:
        name=fields[header.index('Event')]
# an optional list can be used to select specific targets within year
        if targetSelection=='All' or name in targetSelection:
          event={}

          RA=fields[header.index('RA(J2000)')]
          RA=string.split(RA,':')
          RA=15.*(float(RA[0]) + (float(RA[1]) + (float(RA[2])/60.))/60.)
          Dec=fields[header.index('Dec(J2000)')]
          Dec=string.split(Dec,':')
          Dec=np.sign(float(Dec[0])) * (
            abs(float(Dec[0])) + (float(Dec[1])+(float(Dec[2])/60.))/60.)

          event['RA']=RA
          event['Dec']=Dec

          event['fb']=float(fields[header.index('fbl')])
          event['Ftot']=float(fields[header.index('I_bl')])
          event['Fs']=float(fields[header.index('I0')])

          event['t0']=float(fields[header.index('Tmax(HJD)')])
# oof, some trouble here from 2016-45 which has '-' for tE and u0
          try:
            event['tE']=float(fields[header.index('tau')])
          except:
            print('TROUBLE: bad data for tE for',name, \
              fields[header.index('tau')])
            event['tE']=100.
          try:
            event['u0']=float(fields[header.index('umin')])
          except:
            print('TROUBLE: bad data for u0 for',name, \
              fields[header.index('umin')])
            event['u0']=1.
              
          if event['u0']==0:
#          print 'NOTE: u0 is zero.  using 1/A instead'
            event['u0']=1./float(fields[header.index('Amax')])
#          print 'name:',event['name']

          try:
            event['A']=float(fields[header.index('Amax')])
          except:
            event['A']=-123.

          events[name]=event
                  
  print('# of events for',year+':',len(events.keys()))

  return events
#______________________________________________________________________________

def getAllOGLEheaders(options):
  '''
  Read in just the header info for all OGLE alerts: 2015 and onward
  '''
  import os

  events={}
  
  #for year in ['2015','2016','2017','2018','2019','2020']:
  for year in ['2015','2016','2017','2018']:
    
   dir=options['dataDir']+'/ogle/'+year+'/'
   # open the list of all the events and their OGLE single-lens fit params
   filename='lenses.par'

   if not os.access(dir+filename, os.R_OK):
    print('SKIP: no OGLE data for this year yet:',year)
   else:
    file=open(dir+filename,'r')

    for line in file:
      if line.strip()=='':
        #print 'skipping a blank line'
        pass
      else:
        fields=line.split()
        if fields[0]=='Event':
          header=fields
        else:
          name=fields[header.index('Event')]

          event={}

          RA=fields[header.index('RA(J2000)')]
          RA=string.split(RA,':')
          RA=15.*(float(RA[0]) + (float(RA[1]) + (float(RA[2])/60.))/60.)
          Dec=fields[header.index('Dec(J2000)')]
          Dec=string.split(Dec,':')
          Dec=np.sign(float(Dec[0])) * (
            abs(float(Dec[0])) + (float(Dec[1])+(float(Dec[2])/60.))/60.)

          event['RA']=RA
          event['Dec']=Dec

          event['fb']=float(fields[header.index('fbl')])
          event['Ftot']=float(fields[header.index('I_bl')])
          event['Fs']=float(fields[header.index('I0')])

          event['t0']=float(fields[header.index('Tmax(HJD)')])
# oof, some trouble here from 2016-45 which has '-' for tE and u0
          try:
            event['tE']=float(fields[header.index('tau')])
          except:
            print('TROUBLE: bad data for tE for',name, \
              fields[header.index('tau')])
            event['tE']=100.
          try:
            event['u0']=float(fields[header.index('umin')])
          except:
            print('TROUBLE: bad data for u0 for',name, \
              fields[header.index('umin')])
            event['u0']=1.
              
          if event['u0']==0:
#          print 'NOTE: u0 is zero.  using 1/A instead'
            event['u0']=1./float(fields[header.index('Amax')])
#          print 'name:',event['name']

          try:
            event['A']=float(fields[header.index('Amax')])
          except:
            event['A']=-123.

          events[name]=event
                  
    print('# of events up to ',year+':',len(events.keys()))

  return events
#______________________________________________________________________________

def getMOAdata(year,options,targetSelection='All',
               verbose=False):

  import os

  if targetSelection==[]:
    targetSelection='All'
  
#  dir='./data/moa/'+year+'/'
  dir=options['dataDir']+'/moa/'+year+'/'

# load in the index of MOA events
  filename='index'+year+'.dat'
  catalog = open(dir+filename, "r")
  lines = catalog.readlines()
  catalog.close()

  events={}
  for line in lines:
    cols=string.split(line)
    moaname=cols[0]

# now load in the data!
#  for name in targets:
    if targetSelection=='All' or moaname in targetSelection:
      event={'data':{}} 

      event['t0']=float(cols[4])
      event['tE']=float(cols[5])
      event['u0']=float(cols[6])
      event['Ftot']=float(cols[7])  
      event['fb']=1.
      try:
        RA=float(cols[2])
        Dec=float(cols[3])
      except:
#        print 'christ. for some reason this RA/Dec is hr:min:sec', \
#              moashortname,cols[2],cols[3]
        RA=string.split(cols[2],':')
        RA=15.*(float(RA[0]) + (float(RA[1]) + (float(RA[2])/60.))/60.)
        Dec=string.split(cols[3],':')
        Dec=np.sign(float(Dec[0])) * (
          abs(float(Dec[0])) + (float(Dec[1])+(float(Dec[2])/60.))/60.)
#        print '  ',RA,Dec
      event['RA']=RA
      event['Dec']=Dec

      if 'BLG' in moaname:
        moashortname='mb'+moaname[2:4]+moaname[-3:]
      elif 'LMC' in moaname:
        moashortname='ml'+moaname[2:4]+moaname[-3:]
      else:
        exit('ERROR: undefined format for name in index.dat')
#      print moashortname
      filename=moashortname+'.phot'
      print('opening photometry file:',filename)
      
# if the MOA data doesn't exist, grab it from the MOA website
      if not os.access(dir+filename, os.R_OK):
        exit('ERROR: no MOA data for this target.  Run downloadData.py')
      catalog = open(dir+filename, "r")
      lines = catalog.readlines()
      catalog.close()
      photom={'time':[],
              'mag':[],'magerr':[],
              'flux':[],'fluxerr':[]}

      ndrop=0
      for line in lines:
        if line[0]!='#' and line[0]!='<' and string.strip(line):
          cols=string.split(line)
          #print 'line',line
          #print 'cols',cols
          time=float(cols[0]) - options['JDadjust']
          dflux=float(cols[1])
          dfluxerr=float(cols[2])
#          print time,len(photom['time'])
          
# point 1288 of mb15003.phot has zero for the JD time
#          if float(cols[0])==0:
# point 15677 of mb15003.phot is weirder.  -28000 days offset
          if time<-5000.:
            print(' NOTE: dropping a strange photom point',filename,time)
            ndrop+=1
# try dropping points with larger error bars
# hmm, this cuts out everything for mb15310. relax or omit
#          elif dfluxerr>700.:
#          elif dfluxerr>7000.:
#          elif dfluxerr>70000.:
#            print ' NOTE: dropping a poor photom point',filename,time,dflux,dfluxerr
#            ndrop+=1
          else:
            if dfluxerr<=0. or np.isnan(dfluxerr):
              print(' skipping bad MOA errorbar:', \
                    moashortname,time,dfluxerr)
            else:
              photom['time'].append(time)
              photom['flux'].append(dflux)
              photom['fluxerr'].append(dfluxerr)
#            photom['mag'].append(dflux)
#            photom['magerr'].append(dfluxerr)
      photom['time']=np.array(photom['time'])
      photom['flux']=np.array(photom['flux'])
      photom['fluxerr']=np.array(photom['fluxerr'])

      print('# of datapoints, # skipped',len(photom['time']),ndrop)
      #print events
      #print moaname
      events[moaname]=event
      events[moaname]['data']['moa']=photom
      events[moaname]['moaname']=moashortname
    else:
      print('not opening any photometry for this guy',moaname)
      
  return events
#__________________________________________________

def addSpitzer(events,year,options,targets='All'):

  if targets=='All' or targets==[]:
    targets=events.keys()

#  dir='./data/spitzer/'
  dir=options['dataDir']+'/spitzer/'
  
  for name in targets:
    obname='ob'+year[-2:]+name[-4:]
    filename=obname+'.dat'     

    print('Reading file: '+dir+filename+'  ...')
    catalog = open(dir+filename, "r")
    lines = catalog.readlines()
    catalog.close()
    photom={'time':[],'mag':[],'magerr':[],
             'flux':[],'fluxerr':[]}
    for line in lines:
      cols=string.split(line)
      photom['time'].append(float(cols[0])+
                            2450000.-options['JDadjust'])

    mag=float(cols[3])
    err=float(cols[4])
# add flux info here, rather than calc many times
    zeropoint=1
    flux=zeropoint*10.**(-mag/2.5)
# use this factor to convert a mag err to flux err
    errorConvertFactor=0.92
    fluxerr=flux*err*errorConvertFactor

# error is just a photon noise type error itk
# include a systematic error in some fashion
#    magerrmin=0.005
#    err=sqrt(magerrmin**2 + err**2)

    photom['mag'].append(mag)
    photom['magerr'].append(err)
    photom['flux'].append(flux)
    photom['fluxerr'].append(fluxerr)
#    print 'mag err',mag,err

    events[name]['data']['spitzer']=photom

  return events
#__________________________________________________

def getMOAheader(year,options,targetSelection='All'):
  '''
  Read in just the header info for MOA alerts
  The year to be read in is specified.
  Optionally select a batch of targets; otherwise all targets considered.
  '''
  import os
  from math import log10

  if targetSelection==[]:
    targetSelection='All'
  
  dir=options['dataDir']+'/moa/'+year+'/'
  filename='index'+year+'.dat'
  print('opening MOA header:',dir+filename)
  if not os.access(dir+filename, os.R_OK):
    exit('ERROR: no MOA data for this year.  Run downloadData.py')
  file=open(dir+filename,'r')

  events={} 
  for line in file:
    if line.strip()=='':
#      print 'skipping a blank line'
      pass
    else:
      fields=line.split()
      name=fields[0]
      if targetSelection=='All' or name in targetSelection:
        #print fields
        event={}
        try:
          event['RA']=float(fields[2])
          event['Dec']=float(fields[3])
        except:
          print('STUPID: one of the MOA index lines has h:m:s RA/Dec!')
          print('  ',fields[2],fields[3])
          RA = 15.*( float(fields[2][:2]) + \
                     float(fields[2][3:5])/60. + \
                     float(fields[2][6:11])/3600. )
          Dec = float(fields[3][1:3]) + \
                float(fields[3][4:6])/60. + \
                float(fields[3][7:12])/3600.
          if fields[3][0]=='-':
            Dec*=-1
          else:
            exit('ERROR: huh? its not in the south?')
          print('  ',RA,Dec)
          event['RA']=float(RA)
          event['Dec']=float(Dec)
        event['t0']=float(fields[4])
        event['tE']=float(fields[5])
        event['u0']=float(fields[6])
        event['Fs']=float(fields[7])  
        #event['Fb']=float(fields[8])  # ???

        #print 't0,tE',event['t0'],event['tE']
        #print 'u0',event['u0']
        
        events[name]=event
                  
  print('# of events for',year+':',len(events.keys()))

  return events
#______________________________________________________________________________

def getAllMOAheaders(options):
  '''
  Read in the header info for all MOA alerts, 2015-
  Wait actually 2015 header file seems to be missing. no 2015 for now
  '''
  import os
  from math import log10

  events={}
  
  #for year in ['2015','2016','2017','2018','2019','2020']:
  for year in ['2015','2016','2017','2018','2019']:
  #for year in ['2015','2016','2017','2018']:
  #for year in ['2016','2017','2018']:

    dir=options['dataDir']+'/moa/'+year+'/'
    filename='index'+year+'.dat'
    print('opening MOA header:',dir+filename)
    file=open(dir+filename,'r')

    for line in file:
      if line.strip()=='':
        #print 'skipping a blank line'
        pass
      else:
        fields=line.split()
        name=fields[0]

        event={}
        try:
          event['RA']=float(fields[2])
          event['Dec']=float(fields[3])
          print('RA Dec is o.k.',year,fields[1],fields[2],fields[3])
        except:
          print('STUPID: one of the MOA index lines has h:m:s RA/Dec!')
          print('  fields',fields)
          print('  # of fields',len(fields))
          print('  OK the problem is that 2017 index file is suddently messed up.  use the old file.  (but then youre missing 73 newish lightcurves!)')
          print('  ',fields[2],fields[3])
          RA = 15.*( float(fields[2][:2]) + \
                     float(fields[2][3:5])/60. + \
                     float(fields[2][6:11])/3600. )
          Dec = float(fields[3][1:3]) + \
                float(fields[3][4:6])/60. + \
                float(fields[3][7:12])/3600.
          if fields[3][0]=='-':
            Dec*=-1
          else:
            exit('ERROR: huh? its not in the south?')
          print('  ',RA,Dec)
          event['RA']=float(RA)
          event['Dec']=float(Dec)
        event['t0']=float(fields[4])
        event['tE']=float(fields[5])
        event['u0']=float(fields[6])
        event['Fs']=float(fields[7])  
        #event['Fb']=float(fields[8])  # ???

        #print 't0,tE',event['t0'],event['tE']
        #print 'u0',event['u0']
        
        events[name]=event
                  
    print('# of events for',year+':',len(events.keys()))

  return events
#______________________________________________________________________________

def addMOA(events,year,options,targets='All'):
#  execfile('hjd.py')
#  from hjd import helio_jd

  if targets=='All' or targets==[]:
    targets=events.keys()

#  dir='./data/moa/'+year+'/'
  dir=options['dataDir']+'/moa/'+year+'/'

# get the list of OGLE/MOA cross-matches
  filename='moa2ogle_'+year+'.txt'
  if not os.access(dir+filename, os.R_OK):
    exit('ERROR: no MOA data for this year.  Run downloadData.py')
  catalog = open(dir+filename, "r")
  lines = catalog.readlines()
  catalog.close()
  nameConversion={'ogle':[],'moa':[]}    
  for line in lines:
    if line[0]!='#':
      cols=string.split(line)
      if cols[0][:4]=='OGLE':
        oglename=cols[0]
        moaname=cols[2]
      elif cols[0][:3]=='MOA':
        moaname=cols[0]
        oglename=cols[2]
      else:
        print('line:',line)
        exit('ERROR: unexpected format for moa2ogle')
      moashortname='mb'+year[-2:]+moaname[-3:]
#      print moashortname
      nameConversion['ogle'].append(oglename)
      nameConversion['moa'].append(moashortname)

# also load in the index.dat file, to get RA/Dec for HJD conversion
  filename='index'+year+'.dat'
  catalog = open(dir+filename, "r")
  lines = catalog.readlines()
  catalog.close()
  raDecInfo={'moashortname':[],'RA':[],'Dec':[]}
  for line in lines:
    cols=string.split(line)
    moaname=cols[0]
    if 'BLG' in moaname:
      moashortname='mb'+moaname[2:4]+moaname[-3:]
    elif 'LMC' in moaname:
      moashortname='ml'+moaname[2:4]+moaname[-3:]
    else:
      exit('ERROR: undefined format for name in index.dat')
#    print moashortname
    raDecInfo['moashortname'].append(moashortname)
    try:
      RA=float(cols[2])
      Dec=float(cols[3])
    except:
#      print 'christ. for some reason this RA/Dec is hr:min:sec', \
#            moashortname,cols[2],cols[3]
      RA=string.split(cols[2],':')
      RA=15.*(float(RA[0]) + (float(RA[1]) + (float(RA[2])/60.))/60.)
      Dec=string.split(cols[3],':')
      Dec=np.sign(float(Dec[0])) * (
        abs(float(Dec[0])) + (float(Dec[1])+(float(Dec[2])/60.))/60.)
#      print '  ',RA,Dec
    raDecInfo['RA'].append(RA)
    raDecInfo['Dec'].append(Dec)
#  exit()
  
# now load in all the data!
  for name in targets:
    if not 'OGLE-'+name in nameConversion['ogle']:
      print('NOTE: no MOA data for this OGLE event:',name)
#      moashortname=''
#      photom=''
      events[name]['moaname']=''
    else:
      index=nameConversion['ogle'].index('OGLE-'+name)
      moashortname=nameConversion['moa'][index]

      filename=moashortname+'.phot'

# if the MOA data doesn't exist, grab it from the MOA website
      if not os.access(dir+filename, os.R_OK):
        exit('ERROR: no MOA data for this target.  Run downloadData.py')
      catalog = open(dir+filename, "r")
      lines = catalog.readlines()
      catalog.close()
      photom={'time':[],
              'mag':[],'magerr':[],
              'flux':[],'fluxerr':[]}

      ndrop=0
      for line in lines:
        if line[0]!='#':
          cols=string.split(line)
          time=float(cols[0]) - options['JDadjust']
          dflux=float(cols[1])
          dfluxerr=float(cols[2])
#          print time,len(photom['time'])
          
# point 1288 of mb15003.phot has zero for the JD time
#          if float(cols[0])==0:
# point 15677 of mb15003.phot is weirder.  -28000 days offset
          if time<-5000.:
#            print 'NOTE: dropping a strange photom point',filename,time
            ndrop+=1
# try dropping points with larger error bars
# hmm, this cuts out everything for mb15310. relax or omit
#          elif dfluxerr>700.:
          elif dfluxerr>7000.:
#            print 'NOTE: dropping a poor photom point',filename,time
            ndrop+=1
          else:

# convert from JD to HJD.  MJD?
# skip this for now.  it's slow
            if 0:
#              print 'time before conversion',time
              ii=raDecInfo['moashortname'].index(moashortname)
              RA=raDecInfo['RA'][ii]
              Dec=raDecInfo['Dec'][ii]
##              time=helio_jd(time, RA, Dec)
#              print 'ra dec',RA,Dec
#              print 'time after conversion',time
#              print
            if dfluxerr<=0. or np.isnan(dfluxerr):
              print(' skipping bad MOA errorbar:', \
                    moashortname,time,dfluxerr)
            else:
              photom['time'].append(time)
              photom['flux'].append(dflux)
              photom['fluxerr'].append(dfluxerr)
#            photom['mag'].append(dflux)
#            photom['magerr'].append(dfluxerr)
      photom['time']=np.array(photom['time'])
      photom['flux']=np.array(photom['flux'])
      photom['fluxerr']=np.array(photom['fluxerr'])

      print('# of datapoints, # skipped',len(photom['time']),ndrop)
      events[name]['data']['moa']=photom
      events[name]['moaname']=moashortname
      
  return events
#__________________________________________________

def readErrorScaling(eventNames):

  errorScaling={}

  dir='./'
  filename='errorScaling.out'
  if os.access(filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()

    for line in lines:
      cols=string.split(line,'|')
      event=string.strip(cols[0])
      if event in eventNames:
        errorScaling[event]=float(cols[2])
# sometimes failed fits are saved as NaN scaling (2015-BLG-1526)
        if np.isnan(errorScaling[event]):
          errorScaling[event]=1.
      else:
        errorScaling[event]=1.

  else:
    print('NOTE: errorScaling.out file missing?')

  return errorScaling
#__________________________________________________

def startWithPastResults(events,options,filename='results.out'):

  dir='./'
  if os.access(filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()

    for line in lines:
      cols=string.split(line,'|')
      name=string.strip(cols[0])      
#      print name
      
# check that its a full line
#  (sometimes there's cutoff lines in the middle of a run, etc)
      if name in events.keys():
       if len(cols)<9:
        print('skipping a garbage line in results.out ',name,len(cols))
       else:
#        print 'old event values',events[name]['u0'], \
#              events[name]['t0'],events[name]['tE'], \
#              events[name]['Ftot'],events[name]['fb'], \
#              events[name]['Fs']
        try:
#          t0=float(string.split(cols[3],'+-')[0]) \
#              + options['JDadjust']
          t0=float(string.split(cols[3],'+-')[0])

#ASDF: this adjustment is problematic.  redo from stratch
# (bear in mind that OGLE data is full JD, UKIRT is already adjusted)
#          if not options['justUKIRT']:
#            t0+=options['JDadjust']

        except:
          t0=666.
#        print 't0',t0,options['JDadjust']
#        print 't0 adjusted',t0-options['JDadjust']
        
          
# only start with the old t0 if it is near the data range
#  otherwise it is presumably junk that won't help in fitting
# maybe do the same with tE and u0
        avEpoch=np.median(events[name]['data']['ukirt']['time'])
#        avEpoch=np.median(events[name]['data']['ogle']['time']) \
#                 - options['JDadjust']               
#        if (t0-options['JDadjust']>7180. 
#            and t0-options['JDadjust']<7220) \
#            or (t0-options['JDadjust']>1180. 
#                and t0-options['JDadjust']<1220):
        if abs(t0-avEpoch)<1000.:
#          print 'using previous-fit info',name
          events[name]['t0']=t0          
          try:
            events[name]['u0']=float(string.split(cols[1],'+-')[0])
            events[name]['tE']=float(string.split(cols[2],'+-')[0])
          except:
            events[name]['u0']=666.
            events[name]['tE']=666.
        else:
#          print 'using previous baseline info only'
          events[name]['tE']=10.
          events[name]['u0']=1.
          events[name]['t0']= \
               np.median(events[name]['data']['ukirt']['time'])
#               np.median(events[name]['data']['ukirt']['time']) + \
#               options['JDadjust']
        try:
          events[name]['Ftot']=float(string.split(cols[4],'+-')[0])
          events[name]['fb']=float(string.split(cols[5],'+-')[0])
        except:
          events[name]['Ftot']=666.
          events[name]['fb']=1.
        events[name]['Fs']=events[name]['Ftot'] - \
                            2.5*math.log10(events[name]['fb'])
#        print 'new event values',events[name]['u0'], \
#              events[name]['t0'],events[name]['tE'], \
#              events[name]['Ftot'],events[name]['fb'], \
#              events[name]['Fs']

# save the previous-fit chi2 value, as indicator of fit quality
        events[name]['chi2best']=float(cols[-8])
        events[name]['chi2base']=float(cols[-7])
        events[name]['chi2drop1']=float(cols[-6])
        events[name]['npoints']=float(cols[-5])
#        events[name]['chi2med']=float(cols[-4])
#        events[name]['chi2red']=float(cols[-3])
        events[name]['chi2red']=float(cols[-4])

#        print 'at the end',events[name]
#        print
  else:
    print('NOTE: results.out file missing?')

  return events
#__________________________________________________

def loadFitResults(filename,dir='./'):

  print('reading results file:',dir+filename)

  events={}
  if os.access(dir+filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()
    
    for line in lines:
      cols=string.split(line,'|')

      if string.strip(cols[0])=='target':
        print(' skipping the header line of',filename)

      else:
        name=string.strip(cols[0])      
        events[name]={}

        events[name]['u0']=float(string.split(cols[1],'+-')[0])
        events[name]['tE']=float(string.split(cols[2],'+-')[0])
        events[name]['t0']=float(string.split(cols[3],'+-')[0])
        events[name]['Ftot']=float(string.split(cols[4],'+-')[0])
        events[name]['fb']=float(string.split(cols[5],'+-')[0])

        events[name]['u0err']=float(string.split(cols[1],'+-')[1])
        events[name]['tEerr']=float(string.split(cols[2],'+-')[1])
        events[name]['t0err']=float(string.split(cols[3],'+-')[1])
        events[name]['Ftoterr']=float(string.split(cols[4],'+-')[1])
        events[name]['fberr']=float(string.split(cols[5],'+-')[1])

        events[name]['chi2best']=float(cols[6])
        events[name]['chi2base']=float(cols[7])
        events[name]['chi2drop1']=float(cols[8])
        events[name]['npoints']=float(cols[9])
        events[name]['chi2red']=float(cols[10])

  else:
    print(' missing results file',filename)
    exit('ERROR: probably a path/dir problem')
    
  return events
#__________________________________________________

def loadSinFitResults(filename,dir='./'):
# this is the same as loadFitResults, except sin fit instead of mcmc fit

  print('reading results file:',dir+filename)
  
  events={}
  if os.access(dir+filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()
    
    for line in lines:
      cols=string.split(line,'|')

      if string.strip(cols[0])=='target':
        print(' skipping the header line of',filename)

      else:
        name=string.strip(cols[0])      
        events[name]={}

        events[name]['P']=float(string.split(cols[1],'+-')[0])
        events[name]['A']=float(string.split(cols[2],'+-')[0])
        events[name]['F']=float(string.split(cols[3],'+-')[0])
        events[name]['T']=float(string.split(cols[4],'+-')[0])

        events[name]['Perr']=float(string.split(cols[1],'+-')[1])
        events[name]['Aerr']=float(string.split(cols[2],'+-')[1])
        events[name]['Ferr']=float(string.split(cols[3],'+-')[1])
        events[name]['Terr']=float(string.split(cols[4],'+-')[1])

        events[name]['chi2best']=float(cols[5])
        events[name]['chi2base']=float(cols[6])
        events[name]['chi2drop1']=float(cols[7])
        events[name]['npoints']=float(cols[8])
        events[name]['chi2red']=float(cols[9])

  else:
    print(' missing results file',filename)
    exit('ERROR: probably a path/dir problem')
    
  return events
#__________________________________________________

def loadSlopeFitResults(filename,dir='./'):
# this is the same as loadFitResults, except slope fit instead of mcmc fit

  print('reading results file:',dir+filename)
  
  events={}
  if os.access(dir+filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()
    
    for line in lines:
      cols=string.split(line,'|')

      if string.strip(cols[0])=='target':
        print(' skipping the header line of',filename)

      else:
        name=string.strip(cols[0])      
        events[name]={}

        events[name]['m']=float(string.split(cols[1],'+-')[0])
        events[name]['b']=float(string.split(cols[2],'+-')[0])

        events[name]['merr']=float(string.split(cols[1],'+-')[1])
        events[name]['berr']=float(string.split(cols[2],'+-')[1])

        events[name]['chi2best']=float(cols[3])
        events[name]['chi2base']=float(cols[4])
        events[name]['chi2drop1']=float(cols[5])
        events[name]['npoints']=float(cols[6])
        events[name]['chi2red']=float(cols[7])

  else:
    print(' missing results file',filename)
    exit('ERROR: probably a path/dir problem')
    
  return events
#__________________________________________________

def loadGridResults(field,ccd=None,dir='gridresults/',
                    filename='gridresults',
                    fakename=None, Verbose=True):
  from math import log10
  
  events={}

  if filename=='gridresults' or filename=='gridresultsOverlap':
    if options['photomType']=='PSF':
      filename+='PSF'
    filename+=year+'_'+field
    if ccd:
      filename+='_'+ccd
    filename+='.out'
  
  if Verbose:
    print('reading gridresults file:',dir+filename)
  if fakename:
    dir='./'
    filename=fakename
  
  if os.access(dir+filename, os.R_OK):
    catalog = open(dir+filename, 'r')
    lines = catalog.readlines()
    catalog.close()

    for line in lines:
      cols=string.split(line,'|')

      if string.strip(cols[0])=='target':
        #print ' skipping the header line of',filename
        pass
      
# check that its a full line
#  (sometimes there's cutoff lines in the middle of a run, etc)
      elif len(cols)!=12:
        name=string.strip(cols[0])      
        print(' skipping a garbage line ',name,len(cols))

      else:
        name=string.strip(cols[0])      
        events[name]={}
# skip this repeat check for now. slow
#        if not name in events.keys():
#          events[name]={}
##          print 'new guy',name
#        else:
#          print ' repeated guy',name

#        print cols
        events[name]['u0']=float(cols[1])
        events[name]['tE']=float(cols[2])
        events[name]['t0']=float(cols[3])
        events[name]['Ftot']=float(cols[4])
        events[name]['fb']=float(cols[5])
# skip this. slow.
#        events[name]['Fs']=events[name]['Ftot'] - \
#                            2.5*log10(events[name]['fb'])

        events[name]['chi2best']=float(cols[6])
        events[name]['chi2base']=float(cols[7])
        events[name]['chi2drop1']=float(cols[8])
        events[name]['npoints']=float(cols[9])
        events[name]['chi2red']=float(cols[10])

  else:
    print(' no gridresults file for field,ccd=',field,ccd)
    
  return events
#__________________________________________________

def setInitialParams(starname,event,options):
  from math import log10

  params={}
  params['u0']=event['u0']
  params['tE']=event['tE']
#  params['t0']=event['t0'] - options['JDadjust']
  params['t0']=event['t0']
  params['Ftot']=event['Ftot']
  params['fb']=event['fb']
  params['Fs']=event['Fs']
  Fbflux=10.**(-params['Ftot']/2.5) - 10.**(-params['Fs']/2.5)
  if Fbflux<0:
    exit('ERROR: source flux greater than total flux')
  fluxmin=1.e-20
  params['Fb']=-2.5*log10(Fbflux + fluxmin)
#  print 'Fb',params['Fb']
  if params['fb']==1.:
    params['logFrat']=2.
  else:
    params['logFrat']=log10(params['fb']/(1.-params['fb']))
  params['chi2']=1.

  # arbitrary guess at the MOA fit parameters
  #  1/10/19  commented out.  probably not used anymore
  #params['FbMOA']=0.
  #params['FsMOA']=0.

  return params
#_________________________________________________

def getOGLEalertTime(events,options):
  '''
  Read in the OGLE alert times (from Radek via Yossi).
   Years 2011 through 2015 are all included in the single file.
  Cross-match the alert list with the events considered now.
  '''
  
#  dir='data/ogle/'
  dir=options['dataDir']+'/ogle/'

  filename='ews_last_obs_b4_ann.dat'
  file=open(dir+filename,'r')
  
  alerts={'event':[],'time':[]} 
  for line in file:
    fields=line.split()
    alerts['event'].append(fields[0])
    alerts['time'].append(float(fields[1]) + 2450000.)

# find the alert for each event in events
  for event in events:
    if 'OGLE-'+event in alerts['event']:
      ii=alerts['event'].index('OGLE-'+event)
      alertTime=alerts['time'][ii]
    else:
      print('ERROR: this event is not in the alert list',event)
#    events[event]['ogleAlert']=alertTime \
    events[event]['data']['ogle']['alertTime']=alertTime \
                                - options['JDadjust']

  return events
#______________________________________________________________________________

def readOGLEcrosscheck(options,dir='ogleCrossCheck/',
                       onlySameYear=False):
  import string
  
  filename='ogleCrossCheck_'+options['photomType']+options['year']+'.txt'
  file=open(dir+filename,'r')

  events={'ogle':[],'ogleRA':[],'ogleDec':[],
          'ukirt':[],'ukirtRA':[],'ukirtDec':[],
          'ukirtDetect':[], 
          'correctT0':[], 
          #'t0':[],'u0':[],'tE':[],'H':[],'fb':[],
          'offset':[]} 
# uh oh, how to handle the extra info for UKIRT-only guys
#  for now (temporary) they're listed in the .2016 Xmatch file
#  but this messes up 2015, so edit that one too I guess

  for line in file:
    #print 'looping',len(events['ogle']),len(events['ukirt'])
    fields=line.split('|')
    if string.strip(fields[0])=='name':
      header=fields
    else:
# NOTE: usual method of zipping header and data doesnt work
#  since there are two 'name's and two 'RA','Dec's
     if options['year']+'-BLG' in fields[0] or not onlySameYear:
      events['ogle'].append(string.strip(fields[0]))
      if string.strip(header[3])=='name':
        events['ukirt'].append(string.strip(fields[3]))
      elif string.strip(header[1])=='name':
        events['ukirt'].append(string.strip(fields[1]))
      else:
        exit('what column for UKIRT name?')

      if string.strip(header[6])=='offset':
#        events['offset'].append(string.strip(fields[6]))
        events['offset'].append(float(fields[6]))
      else:
        exit('what column for offset?')
        
      events['ogleRA'].append(float(fields[1]))
      events['ogleDec'].append(float(fields[2]))
      events['ukirtRA'].append(float(fields[4]))
      events['ukirtDec'].append(float(fields[5]))

# this is crappy. will be cleaned up later...
# 9/5/17 actually, just remove the whole thing. it's not used at all
#      try:
#        events['t0'].append(float(fields[8]))
#        events['u0'].append(float(fields[9]))
#        events['tE'].append(float(fields[10]))
#        events['H'].append(float(fields[11]))
#        events['fb'].append(float(fields[12]))
#      except:
#        events['t0'].append(fields[8])
#        events['u0'].append(fields[9])
#        events['tE'].append(fields[10])
#        events['H'].append(fields[11])
#        events['fb'].append(fields[12])

      if string.strip(header[7])=='ukirtDetect':
        try:
          events['ukirtDetect'].append(float(fields[7]))
        except:
          events['ukirtDetect'].append(fields[7])
      else:
        exit('update UKIRT_OGLE_Xmatch to include UKIRT detection')

      if string.strip(header[8])=='correctT0':
        try:
          events['correctT0'].append(float(fields[8]))
        except:
          events['correctT0'].append(fields[8])
      else:
        print(string.strip(header[8]))
        exit('update UKIRT_OGLE_Xmatch to include T0 check')

  return events
#______________________________________________________________________________

def readMOAcrosscheck(options, onlySameYear=False):

  # THIS IS IDENTICAL TO readOGLEcrosscheck() EXCEPT FOR A FEW 'ogle --> moa'

  import string

  dir='moaCrossCheck/'
  filename='moaCrossCheck_'+options['photomType']+options['year']+'.txt'
  file=open(dir+filename,'r')

  events={'moa':[],'moaRA':[],'moaDec':[],
          'ukirt':[],'ukirtRA':[],'ukirtDec':[],
          'ukirtDetect':[], 
          'correctT0':[], 
          #'t0':[],'u0':[],'tE':[],'H':[],'fb':[],
          'offset':[]} 
# uh oh, how to handle the extra info for UKIRT-only guys
#  for now (temporary) they're listed in the .2016 Xmatch file
#  but this messes up 2015, so edit that one too I guess

  for line in file:
    #print 'looping',len(events['moa']),len(events['ukirt'])
    fields=line.split('|')
    if string.strip(fields[0])=='name':
      header=fields
    else:
     if options['year'] in fields[0] or not onlySameYear:
# NOTE: usual method of zipping header and data doesnt work
#  since there are two 'name's and two 'RA','Dec's
      events['moa'].append(string.strip(fields[0]))
      if string.strip(header[3])=='name':
        events['ukirt'].append(string.strip(fields[3]))
      elif string.strip(header[1])=='name':
        events['ukirt'].append(string.strip(fields[1]))
      else:
        exit('what column for UKIRT name?')

      if string.strip(header[6])=='offset':
#        events['offset'].append(string.strip(fields[6]))
        events['offset'].append(float(fields[6]))
      else:
        exit('what column for offset?')
        
      events['moaRA'].append(float(fields[1]))
      events['moaDec'].append(float(fields[2]))
      events['ukirtRA'].append(float(fields[4]))
      events['ukirtDec'].append(float(fields[5]))

# this is crappy. will be cleaned up later...
#      try:
#        events['t0'].append(float(fields[8]))
#        events['u0'].append(float(fields[9]))
#        events['tE'].append(float(fields[10]))
#        events['H'].append(float(fields[11]))
#        events['fb'].append(float(fields[12]))
#      except:
#        events['t0'].append(fields[8])
#        events['u0'].append(fields[9])
#        events['tE'].append(fields[10])
#        events['H'].append(fields[11])
#        events['fb'].append(fields[12])

      if string.strip(header[7])=='ukirtDetect':
        try:
          events['ukirtDetect'].append(float(fields[7]))
        except:
          events['ukirtDetect'].append(fields[7])
      else:
        exit('update UKIRT_MOA_Xmatch to include UKIRT detection')

      if string.strip(header[8])=='correctT0':
        try:
          events['correctT0'].append(float(fields[8]))
        except:
          events['correctT0'].append(fields[8])
      else:
        exit('update UKIRT_MOA_Xmatch to include T0 check')
        
  return events
#______________________________________________________________________________

def getOGLEvars():
  import string

  dir='./'
  filename='ogleVars.txt'
  file=open(dir+filename,'r')

  varStars={'name':[],
            'type':[],'subtype':[],
            'Vmag':[],'Imag':[],
            'ra':[],'dec':[]} 
  for line in file:
    if line[0]=='#':
      pass
#      print 'skipping header line'
    else:
      fields=line.split()
      varStars['name'].append(string.strip(fields[0]))
      varStars['type'].append(string.strip(fields[5]))
      varStars['subtype'].append(string.strip(fields[6]))
      varStars['Vmag'].append(float(fields[8]))
      varStars['Imag'].append(float(fields[7]))
      ra_hms=string.strip(fields[3])
      ra_hr=float(ra_hms[:2])
      ra_min=float(ra_hms[3:5])
      ra_sec=float(ra_hms[6:])
#      print 'ra',ra_hr,'x',ra_min,'x',ra_sec
      ra=ra_hr+ (ra_min + (ra_sec/60.))/60.
      ra=ra/24.*2.*math.pi
      varStars['ra'].append(ra)

      dec_hms=string.strip(fields[4])
      dec_deg=float(dec_hms[:3])
      dec_min=float(dec_hms[4:6])
      dec_sec=float(dec_hms[7:])
#      print 'dec',dec_deg,'x',dec_min,'x',dec_sec
      dec=abs(dec_deg)+ (dec_min + (dec_sec/60.))/60.
      dec=dec*np.sign(dec_deg)
      dec=dec/360.*2.*math.pi
      varStars['dec'].append(dec)

#      print ra,dec
#  exit()
  varStars['ra']=np.array(varStars['ra'])
  varStars['dec']=np.array(varStars['dec'])
# these strings have to be numpy arrays too, for whereing
  varStars['name']=np.array(varStars['name'])
  varStars['type']=np.array(varStars['type'])
  varStars['subtype']=np.array(varStars['subtype'])
  return varStars
#______________________________________________________________________________


def getVarStarlist():
  import string

  dir='./'
  filename='UKIRT_OGLE_VARmatch.txt.temp'
  file=open(dir+filename,'r')

  varStars={'name':[],'oglename':[],
            'type':[],'subtype':[],
#            'Vmag':[],'Imag':[],
            'oglera':[],'ogledec':[],
            'ukirtra':[],'ukirtdec':[],
            'offset':[]}
  for line in file:
    fields=line.split('|')
    if string.strip(fields[0])=='name':
      pass
#      print 'skipping header line'
    else:
      varStars['oglename'].append(string.strip(fields[0]))
      varStars['type'].append(string.strip(fields[1]))
      varStars['subtype'].append(string.strip(fields[2]))
      varStars['name'].append(string.strip(fields[3]))
# no need to read in RA/Dec info
#      varStars['oglera'].append(string.strip(fields[4]))
#      varStars['ogledec'].append(string.strip(fields[5]))
#      varStars['ukirtra'].append(string.strip(fields[6]))
#      varStars['ukirtdec'].append(string.strip(fields[7]))
      varStars['offset'].append(float(fields[8]))

  return varStars
#______________________________________________________________________________


def readFirstpassUKIRT(year,photomType='CASU',dir='firstpass/',
                       fieldSelect='all',ccdSelect='all',
                       altBand=False,oldData=False):
  import csv,string
  header=''

# special case to load in the previous data reduction, for cross-match IDs 
  if oldData:
    dir += '2017partial/'
    if year!='2017':
      exit('oldData is set up for the 2017 final data download')

  fields,ccds = allFields(year)

  if altBand and (year=='2015' or year=='2016'):
    print('   (this year has no altband)',year)
    return {'RA':[],'Dec':[]}
  
  data={}
  for field in fields:
    if field==fieldSelect  or  fieldSelect=='all':
        filename='firstpass_'+photomType+year+'_'+field
        if altBand:
          filename+='_band2'
        filename+='.txt'
        
        print('opening a firstpass file:',filename)
        file=open(dir+filename,'r')
        file=csv.reader(file,delimiter='|')

        header=next(file)
        for i in range(len(header)):
          header[i]=string.strip(header[i])

        if data=={}:
# only initialize the data dictionary for the first file
          for field in header:    
            data[field]=[]
          firstheader=header
        else:
# make sure that subsequent files have the same header as previous one
#  (these two checks are redundant)
          if header!=firstheader:
            print('TROUBLE: headers dont match',header,firstheader)
            exit('ERROR: headers dont match! 1')
          for field in header:    
            if not field in data.keys():
              print('field:',field,'  keys:',data.keys())
              exit('ERROR: headers dont match! 2')
        
        for line in file:
          ccdHeaderPosition=header.index('ccd')
          ccdNow=string.strip(line[ccdHeaderPosition])
          #print 'ccd',ccdNow

          # special option: only select a particular chip
          if ccdNow==ccdSelect or ccdSelect=='all':
            for field,value in zip(header,line):
              if string.strip(field):
#                print field,'value:',value
# float() doesn't work for the 'name' field
                try:
                  data[field].append(float(value))
                except:
                  data[field].append(string.strip(value))

        print('  number of data now:',len(data['name']))
#      else:
#        print ' skipping this field:',fieldNow

# convert slow lists to fast arrays
  for field in header:
    data[field]=np.array(data[field])

# the final '|' creates a blank field at the end.  remove it
  if '' in data.keys():
    data.pop('')
        
  return data
#______________________________________________________________________________

def readSecondPass(year,photomType='CASU',dir='secondpass/'):
  #                   altBand=False):
  # always read in the altBand, if there is one
  
  import csv,string
  header=''

  filename='secondpass_'+photomType+year+'.txt'
        
  print('opening a secondpass file:',filename)
  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter='|')

  header=next(file)
  for i in range(len(header)):
    header[i]=string.strip(header[i])
  print('  header:',header)
  
  data={}
  for field in header:    
    data[field]=[]
    
  for line in file:      
    for field,value in zip(header,line):
      #print 'field,value',field,value
      if field=='field' or field=='ccd':
        data[field].append(string.strip(value))
      else:
        try:
          data[field].append(float(value))
        except:
          data[field].append(string.strip(value))

  # always read in the altBand, if there is one
  altfilename='secondpass_'+photomType+year+'_band2.txt'
        
  print('opening an altband secondpass file:',filename)
  file=open(dir+altfilename,'r')
  file=csv.reader(file,delimiter='|')

  altheader=next(file)
  for i in range(len(altheader)):
    altheader[i]=string.strip(altheader[i])
  print('  altheader:',altheader)
  if altheader!=header:
    exit('ERR: altband file has a different header')
  
  for line in file:      
    for field,value in zip(header,line):
      #print 'field,value',field,value
      if field=='field' or field=='ccd':
        data[field].append(string.strip(value))
      else:
        try:
          data[field].append(float(value))
        except:
          data[field].append(string.strip(value))


  # convert slow lists to fast arrays
  for field in header:
    data[field]=np.array(data[field])

  # the final '|' creates a blank field at the end.  remove it
  if '' in data.keys():
    data.pop('')
        
  return data
#______________________________________________________________________________

def readOGLEfields(filename,dir='fieldLocations/'):
  import csv,string

  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter=' ')

  header=next(file)
  for i in range(len(header)):
    header[i]=string.strip(header[i])
  
  data={}
  for field in header:    
    data[field]=[]
  for line in file:
    for field,value in zip(header,line):
      if string.strip(field):
#        print field,'value:',value
        try:
          data[field].append(float(value))
        except:
          data[field].append(value)
#    print data
#    print 'line',line
#    exit()

# convert coordinates to RA/Dec degrees
  if 'R.A.' in data:
    dechms=data['Dec']
    data['RA']=[]
    data['Dec']=[]
# awkward - 'Dec' already exists!
    for RAhms,Dechms in zip(data['R.A.'],dechms):
      RA=RAhms
      RA,Dec=convertRADec(RAhms,Dechms)
#      print 'RA:',RA,RAhms
#      print 'Dec:',Dec,Dechms
      data['RA'].append(RA)
      data['Dec'].append(Dec)
  else:
# awkward - 'RA' already exists!
    rahms=data['RA']
    data['RA']=[]
    data['Dec']=[]
    for RAhms,Dechms in zip(rahms,data['DEC']):
      RA=RAhms
      RA,Dec=convertRADec(RAhms,Dechms)
#      print 'RA:',RA,RAhms
#      print 'Dec:',Dec,Dechms
      data['RA'].append(RA)
      data['Dec'].append(Dec)
      
  for field in data:
    data[field]=np.array(data[field])

  return data
#______________________________________________________________________________

def readEpochRange(options):
  import csv,string

  year=options['year']
  photomType=options['photomType']

# pick some representative data file for the overall time range.
#  note that these HJD files do vary a bit from field to field
#  also they vary a lot by band, H vs K
  if year=='2015' or year=='2016':
    filename = photomType+'_hjd_'+year[-2:]+'_11_1.dat'
  else:
    filename = photomType+'_hjd_'+year+'_c1_1_1_K.dat'

  dir = options['dataDir']+'/ukirt/'+year+'/'

  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter=' ')
  
  data=[]
  for line in file:
    data.append(float(line[0]))
  data=np.array(data)
  
  return data
#______________________________________________________________________________

def readEpochbyField(options,field):
  import csv,string
  
  year=options['year']
  photomType=options['photomType']

  if year=='2015' or year=='2016':
    filename = photomType+'_hjd_'+year[-2:]+'_'+field+'_1.dat'
  else:
    if 'c' in field or 's1' in field or 's2' in field:
      filename = photomType+'_hjd_'+year+'_'+field+'_1_K.dat'
    else:
      filename = photomType+'_hjd_'+year+'_'+field+'_1_H.dat'

  dir = options['dataDir']+'/ukirt/'+year+'/'

  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter=' ')  
  
  data=[]
  for line in file:
    data.append(float(line[0]))
  data=np.array(data)
  
  return data
#______________________________________________________________________________

def readFakedata(filename='fakeData.out',dir='./'):
                                                     
  infile=open(dir+filename,'r')
  t0s=pickle.load(infile)
  tEs=pickle.load(infile)
  u0s=pickle.load(infile)
  Nrep=pickle.load(infile)
  epochs=pickle.load(infile)
  fakedati=pickle.load(infile)

# put the fake data into the same format as real data
#  no need to save param info here itk
#  just unwrap fakedati and put it as time/flux
#  for t0 in t0s:
#    for tE in tEs:
#      for u0 in u0s:
  dims=fakedati.shape
  print('dims',dims)
  if len(dims)!=5:
    print('ERROR: strange shape for fakedati',dims)
  mags=np.reshape(fakedati,(dims[0]*dims[1]*dims[2]*dims[3],dims[4]))

  return t0s,tEs,u0s,Nrep,epochs,mags
#______________________________________________________________________________

def convertFakeEvents(epochs,mags):
  events={}
  ii=0
  for mag in mags:
    ii+=1
    eventName='fakedata'+str('%3i' %ii)
#    print eventName
               
    event={'u0':1.,'tE':10.,
           't0':np.median(epochs),
           'Ftot':np.median(mag),
           'Fs':np.median(mag),
           'fb':1,
           'data':{'ukirt':{}}}
    event['data']['ukirt']['time']=epochs
    event['data']['ukirt']['mag']=mag
    event['data']['ukirt']['magerr']=np.array([0.02]*len(epochs))

    errorConvertFactor=0.92
    zeropointK=1.  # this really has to be 1 or there's offset
    fluxes=zeropointK*10.**(-mag/2.5)
    event['data']['ukirt']['flux']=fluxes
    event['data']['ukirt']['fluxerr']= \
                      event['data']['ukirt']['magerr'] \
                      *fluxes*errorConvertFactor
#    print 'fake event',event
#    exit()

    events[eventName]=event
  return events    

#______________________________________________________________________________

def readUKIRTdata(year,field,options,targetSelection='All',ccd=None,
                  altBand=False,specialC4_2=False,
                  useBools=True,includeGlitches=False):
  '''
  Read in the UKIRT photometry from the /data directory.
  The year/field/ccd is specified directly (** do not use the values in 'options' **)
  The type of photometry (CASU vs PSF) is set in 'options'.
  Optionally select a batch of targets; otherwise all targets considered.
  '''
  import os

  #****
  #if options['photomType']=='PSF':
  #  useBools = False
  #useBools = False
  #****

  if options['photomType']=='CASU':
    dir=options['dataDir']+'/ukirt/'+year+'/'
    photomType='CASU'
  elif options['photomType']=='PSF':
    # use the alternate path to data for PSF (it's too big to store locally)
    dir=options['dataDirPSF']+'/ukirt/PSF'+year+'/'
    if options['dataDirPSF']=='/Volumes/Data/':
      dir=options['dataDirPSF']+'/ukirtlightcurves/PSF'+year+'/'
    photomType='PSF'
  else:
    exit('ERROR: unknown photomType in readUKIRTdata')

# specify the waveband
  primaryband,secondaryband = wavebands(year,field)
  if altBand:
    if secondaryband:
      band='_'+secondaryband
    else:
      # there is no alternate band for 2015 and 2016.  return empty
      return {}
  else:
    band='_'+primaryband
  savedboolband = band
  if year=='2015' or year=='2016':
    # 4/23/19 PSF photometry redone for 2015/2016
    #   has same naming scheme as 2017- now, so drop the band='' line
    #   oop. need to keep it for CASU still, so add this conditional:
    if options['photomType']!='PSF':
      band = ''
    shortyear=year[-2:]
  else:
    shortyear = year  # the file names from yossi have all 4 digits in 2017
    # savannah format
    #if options['photomType']=='PSF' and year=='2017':
    if options['photomType']=='PSF':
      shortyear=year[-2:]
  #print '  field,band,shortyear:',field,band,shortyear
  
  # special for savannah format
  #if options['photomType']=='PSF' and year=='2017':
  #  if altBand:
  #    dir=dir+secondaryband+'/'
  #  else:
  #    dir=dir+primaryband+'/'
  #print 'dir',dir
  
# get the list of files in the directory
  filelist=os.listdir(dir)
  #print 'UKIRT files:',filelist
  #print ' # of UKIRT files:',len(filelist)
  #exit('testread')
  
# determine which ccd's have data in this field
  ccdlist=[]
  #print filelist
  for i in range(10):
    #print 'checking:',photomType+'_hjd_'+shortyear+'_'+field+'_'+str(i)+band+'.dat'
    if photomType+'_hjd_'+shortyear+'_'+field+'_'+str(i)+band+'.dat' in filelist:
      #print ' found it',i
      ccdlist.append(str(i))
    #else:
    #  print 'missing:',photomType+'_hjd_'+shortyear+'_'+field+'_'+str(i)+band+'.dat'
  if ccdlist==[]:
    print('  looking for:',photomType+'_hjd_'+shortyear+'_'+field+'_'+str(i)+band+'.dat')
    # northern fields are only H - no K - in 2018
    if altBand and year=='2018' and 'n' in field:
      print('NO altBand data for:',year,field)
    else:
      exit('ERROR: no ccds available for this field (in readData.py)')
  # allow ccd to be specified in the call (if there are any ccds available)
  elif ccd:
    #print 'ccd',ccd,str(ccd)
    ccdlist=[str(ccd)]
  #print ' CCD list for this field:',field,ccdlist
  
# use this factor to convert a mag err to flux err
  errorConvertFactor=0.92
# set the I band zeropoint, to convert magnitudes to Jy
  zeropointI=1.

  events={} 

  iskip=0
  ionly=-666
# asdf:  use this to just do a few lines, for testing
#  ionly=5
#  ionly=10000
    
# skip some, to get past a bunch of 1-epoch crap (particularly for the PSF reduction)
#  iskip=1000
  ionly+=iskip
  if targetSelection!='All':
    ionly=-666
    iskip=0
#  iskip=69000
#  ionly=iskip+2000
#  iskip=28000
#  ionly=iskip+1000
  
  for ccd in ccdlist:
    fileext = shortyear+'_'+field+'_'+ccd+band
    fileextbool = field+'_'+ccd+savedboolband  # 2015/2016 bool still has the band in name 
    
    #print 'OPENING file extension:',fileext
# first get the observing times (for all events) 
    file=open(dir+photomType+'_hjd_'+fileext+'.dat','r')
    lines=file.readlines()
    file.close()
    times=[]
    for line in lines:     
      times.append(float(line))
    times=np.array(times)
    print(' # of observing epochs',len(times),field+'_'+ccd+band)
    
    radecfilename=photomType+'_'+fileext+'.cat'
    file=open(dir+radecfilename,'r')
    radeclines=file.readlines()
    if photomType=='PSF':
      otherinfofilename=photomType+'_'+fileext+'.catadd'
      file=open(dir+otherinfofilename,'r')
      otherinfolines=file.readlines()
    else:
      otherinfolines=['']*len(radeclines)

# then get the fluxes/errors for each event
    file = open(dir+photomType+'_lc_'+fileext+'.dat','r')
    maglines = file.readlines()
    file.close()
    file = open(dir+photomType+'_err_'+fileext+'.dat','r')
    errlines = file.readlines()
    file.close()

    # get the boolean flags for bad epochs for this field+ccd
    #  oops, these are missing the year in the name,
    #  so different fileext definition needed
    #if year=='2018':
    # 2/11/19 should have bool files for all years now
    #if useBools and (year=='2018' or year=='2015' or year=='2016'):
    if useBools:
      file = open(dir+photomType+'_bool_'+fileextbool+'.dat','r')
      goodFlags = file.readlines()
      #print goodFlags
      goodFlags = np.array(goodFlags,dtype=int)
      #print goodFlags
      file.close()
      print(' # of epochs in the bad epoch flag list:',len(goodFlags))
      if len(goodFlags) != len(times):
        exit('ERROR: wrong length for the bool file')
    else:
      # if there's no bool file, then just keep everything of course
      goodFlags = np.array([1]*len(times))

    # 2/11/19 (for calculating the secondpass info, need to keep all frames)
    if includeGlitches:
      goodFlags = np.array([1]*len(times))

    iname=0
    for magline,errline,radecline,otherinfoline in zip(
      maglines,errlines,radeclines,otherinfolines):

      # skip the header row (it's only there for the PSF photometry)
      if magline[:3]=='col':
        # verify that all files have a header row
        if errline[:3]!='col':
          exit('no header for errs?')
        elif radecline[:3]!='ALP':
          exit('no header for radec?')
        elif otherinfoline[:3]!='col':
          exit('no header for otherinfo?')
      else:
        iname+=1
        if iname<=iskip:
          pass
        else:
          # old naming scheme didn't include year or photomType
          #name='UK_'+field+'_'+ccd+'_'+str(iname)
          name='UK'+year+'_'+field+'_'+ccd+band+'_'+photomType[0]+str(iname)
          #print
          #print 'reading in:',name
 
          event={'name':name,
                 'data':{},  
                 'field':field,'ccd':ccd,'band':band}

          radec=string.split(radecline)
          otherinfo=string.split(otherinfoline)
              
          if photomType=='CASU':
            if len(radec)!=7:
              exit('ERROR: should be 7 columns for RA/Dec cat')

            event['ra']=float(radec[1])
            event['dec']=float(radec[2])

            # also save the other columns in the .cat file, for validation
            event['namecheck']=radec[0]
            event['avmagcheck']=radec[3]
            event['dispcheck']=radec[4]
            event['numcheck']=radec[5]
            
          elif photomType=='PSF':
            if len(radec)!=2:
              exit('ERROR: should be 2 columns for RA/Dec cat')
            elif len(otherinfo)!=5:
              exit('ERROR: should be 5 columns for otherinfo catadd')

            event['ra']=float(radec[0])
            event['dec']=float(radec[1])

            # also save the other columns in the .cat/.catadd files, for validation
            event['namecheck']='no name info'
            event['avmagcheck']=otherinfo[0]
            event['dispcheck']=otherinfo[1]
            event['numcheck']=otherinfo[2]

          else:
            exit('ERROR: bad photomType in readUKIRTdata')
                     
          mags=np.array(string.split(magline),dtype=float)
          errors=np.array(string.split(errline),dtype=float)
          if len(mags)!=len(times):
            print(len(mags),len(times))
            print('ERROR: not enough flux entries')
            exit('fix this1!')
          elif len(errors)!=len(times):
            print('ERROR: not enough flux entries')
            exit('fix this2!')
            
#  ****** SYSTEMATIC ERROR INCLUDED HERE!!! ******            
          systematicErr=0.01
          errors=np.sqrt(errors**2 + systematicErr**2)
           
          #photom={'time':times,
          #        'mag':mags,'magerr':errors}
          #
          # 10/31/18 REMOVE BAD EPOCHS!!
          #print ' # of times len change:',len(times), \
          #  len(times[np.where(goodFlags==1)])
          photom = {'time':times[np.where(goodFlags==1)],
                    'mag':mags[np.where(goodFlags==1)],
                    'magerr':errors[np.where(goodFlags==1)]}
          
# drop NaN crap; causes NaN loglikelihoods
# careful! for the gridFit method, you really want to keep time grid
#
# this is not needed anymore.  equivalent text is put at start of MCMC fitting
# a but wait, will this be redundantly done by sinfit etc?
# ASDF: CHECK THIS LATER...
#          if not options['gouldStyle']:
#            goodEpochs=np.where(np.isfinite(photom['magerr']))
#            photom['time']=photom['time'][goodEpochs]
#            photom['mag']=photom['mag'][goodEpochs]
#            photom['magerr']=photom['magerr'][goodEpochs]

# convert magnitudes to fluxes            
          photom['flux']=zeropointI * 10.**(-photom['mag']/2.5)
          photom['fluxerr']=photom['magerr']* \
                             photom['flux']*errorConvertFactor
            
          event['data']['ukirt']=photom
                  
# set some arbitrary initial values for the fitting routine
#  (not needed for Gould-style, but MCMC has to start somewhere)
# move to inside of the eventfitter region
# capability to restart with old params is now removed itk
#          event['tE']=10.
#          event['u0']=1.
#          event['t0']=np.median(photom['time']) 
#          event['fb']=0.9
#          event['Ftot']=np.median(photom['mag'])
#          event['Fs']=np.median(photom['mag'])

          if specialC4_2 and field=='c4_2' and year=='2017':
            #print field
            #print ccd
            #print event.keys()
            #print event['data']['ukirt'].keys()
            #print 'npoints was',len(event['data']['ukirt']['flux'])
            normalData = np.where(event['data']['ukirt']['time'] < 8000.)
            for key in event['data']['ukirt'].keys():
              #print 'fixing',key
              event['data']['ukirt'][key] = event['data']['ukirt'][key][normalData]
            #print 'npoints is ',len(event['data']['ukirt']['flux'])

            # actually npoints is not defined in event yet
            #print 'npoints was',event['data']['ukirt']['npoints']
            #event['data']['ukirt']['npoints'] = len(event['data']['ukirt']['flux'])
            #print 'npoints is ',event['data']['ukirt']['npoints']

            #exit('special data removal')

          events[name]=event

      if iname==ionly:
        return events

  print(' # of events for',year,field,':',len(events.keys()))
   
# optionally limit the event list to some specific batch of targets
  if targetSelection!='All' and targetSelection!=[]:
    allevents = events
    print('   culling the event list',len(events))
    if len(allevents)<10:
      print(allevents)
    events = {}
    for name in targetSelection:
#      if name:  # special check added for cases where altname is blank
# no this fails the check below.  just don't pass blank ('') in here
        print('targetselection name',name,len(allevents),len(events))
        events[name] = allevents[name]
    print('   culled the event list',len(events))
    if len(events)!=len(targetSelection):
      print('NOTE: not all targets found!',len(targetSelection))
      exit('testing the target culling')

  return events
#______________________________________________________________________________

def readUKIRTbyeye(filename='byEye_RADec.txt', dir='byEye/'):
# this is the improved version, where the file contains RA/Dec info as well as name
  
  file=open(dir+filename,'r')
  
  data={}
  for line in file:
    if line.strip()!='':
      fields=line.split('|')

      if string.strip(fields[0])=='name':
        header=[]
        for field in fields:
          header.append(string.strip(field))
        for column in header:
          data[column]=[]
      else:
        if len(fields)==1:
          #print 'skipping some notes:',string.strip(line)
          for column in header:
            data[column].append('only a note')
          #print data
          data['notes'][-1] = string.strip(fields[0])
        else:
          for column in header:
            if column=='RA' or column=='Dec':
              data[column].append(float(fields[header.index(column)]))
            else:
              data[column].append(string.strip(fields[header.index(column)]))
              
  return data
#______________________________________________________________________________

def readUKIRTbyeyeNAMESONLY(filename='byEyeGuys.txt', dir='byEye/'):
  import csv

  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter=' ')
  
  data=[]
  for line in file:
    data.append(string.strip(line[0]))
  
  return data
#______________________________________________________________________________

def readBADogle(filename='badOGLE.txt', dir='ogleCrossCheck/'):
  import csv

  file=open(dir+filename,'r')
  file=csv.reader(file,delimiter=' ')
  
  data=[]
  for line in file:
    data.append(string.strip(line[0]))
    
  return data
#______________________________________________________________________________

def getGOODogle(options,
                offsetLimit_arcsec=3.,
                dropWrongT0=True):
  
  ogleEvents = readOGLEcrosscheck(options)
  
  Nogle = len(ogleEvents['ogle'])
  print('# of OGLE events in the field',Nogle)

  badOGLE = readBADogle()
       
  goodOGLE = {'ukirtname':[],'oglename':[],
              'ukirtDetect':[],'year':[],
              'RA':[],'Dec':[]}
  for i in range(Nogle):
    if ogleEvents['ogle'][i] in badOGLE:
        #print 'DROP: bad OGLE (duplicate or CV)',ogleEvents['ogle'][i]
        pass
    elif ogleEvents['offset'][i]>offsetLimit_arcsec:
        #print 'DROP: not a real OGLE match',ogleEvents['offset'][i],ogleEvents['ukirt'][i]
        pass
    elif ogleEvents['correctT0'][i]==0 and dropWrongT0:
        #print 'DROP: event not during UKIRT observing',ogleEvents['ogle'][i]
        pass
    else:
        goodOGLE['ukirtname'].append(ogleEvents['ukirt'][i])
        goodOGLE['oglename'].append(ogleEvents['ogle'][i])
        if ogleEvents['ukirt'][i][2:6]!=options['year']:
            exit('ERROR: name doesnt match year')
        goodOGLE['year'].append(options['year'])
        goodOGLE['ukirtDetect'].append(ogleEvents['ukirtDetect'][i])
        #print ogleEvents.keys()
        goodOGLE['RA'].append(ogleEvents['ukirtRA'][i])
        goodOGLE['Dec'].append(ogleEvents['ukirtDec'][i])
  #Nogle = len(goodOGLE['ukirtname'])
  #print '# of GOOD OGLE events with right t0',Nogle

  return goodOGLE
#______________________________________________________________________________

def getGOODmoa(options,
               offsetLimit_arcsec=3.,
               dropWrongT0=True):

  moaEvents = readMOAcrosscheck(options)
  Nmoa = len(moaEvents['moa'])
  print('# of MOA events in the field',Nmoa)
  
  goodMOA = {'ukirtname':[],'moaname':[],
              'ukirtDetect':[],'year':[],
              'RA':[],'Dec':[]}
  for i in range(Nmoa):
    if moaEvents['offset'][i]>offsetLimit_arcsec:
      #print 'DROP: not a real MOA match',moaEvents['offset'][i],moaEvents['ukirt'][i]
      pass
    elif moaEvents['correctT0'][i]==0 and dropWrongT0:
      #print 'DROP: MOA event not during UKIRT observing',moaEvents['moa'][i]
      pass
    else:
      #print 'new MOA time',moaEvents['ukirt'][i]
      goodMOA['ukirtname'].append(moaEvents['ukirt'][i])
      goodMOA['moaname'].append(moaEvents['moa'][i])
      goodMOA['year'].append(options['year'])
      goodMOA['ukirtDetect'].append(moaEvents['ukirtDetect'][i])
      goodMOA['RA'].append(moaEvents['ukirtRA'][i])
      goodMOA['Dec'].append(moaEvents['ukirtDec'][i])
  #Nmoa = len(goodMOA['ukirtname'])
  #print 'final # of events',Nmoa

  return goodMOA
#______________________________________________________________________________

#def readStats(filename='summaryInfo.txt', dir='summaryStats/'):
def readStats(filename='statsTable.tex', dir='summaryStats/'):

# old format:
#  2015 & all         & 4.5 &  6,669,239 & 562 & 14 &  9 & 8108 & TBD & 33 \\
# new format:
#  2015 & all     & 3.4 & 50 & 146 & 6,669,239 & 562 & 13 &  9 & 8108 & TBD & \\

  file=open(dir+filename,'r')
    
  data=[]
  for line in file:
    cols = string.split(line,'&')
    
    #print len(cols)
    if len(cols)!=12:
      #print 'skipping table wrapping row',line
      pass
    else:
      data.append({})
      data[-1]['year']=string.strip(cols[0])
      data[-1]['field']=string.strip(cols[1])
      data[-1]['area']=float(cols[2])
      data[-1]['Nlc']=float(string.replace(cols[5],',',''))
      data[-1]['chi500']=float(cols[6])
      data[-1]['MLclasschi500']=float(cols[7])
      data[-1]['MLclassOGLEchi500']=float(cols[8])    
      data[-1]['chi100']=float(cols[9])
      data[-1]['MLclasschi100']=float(cols[10])
      data[-1]['MLclassOGLEchi100']=float(string.replace(cols[11],'\\',''))
      #data[-1]['byeye']=float(cols[7])
      #data[-1]['byeyeOGLE']=float(cols[8])    


      if data[-1]['field']=='b=$-1.0^\circ$' or \
         data[-1]['field']=='b=$-1.7^\circ$' or \
         data[-1]['field']=='b=$-2.4^\circ$' or \
         data[-1]['field']=='b=$-3.0^\circ$':
        data[-1]['year'] = '2016'
      elif data[-1]['field']=='C1-C4' or \
           data[-1]['field']=='S1-S2' or \
           data[-1]['field']=='S3-S6' or \
           data[-1]['field']=='N1-N4': 
        data[-1]['year'] = '2017'

      if data[-1]['year']=='':
        exit('ERR: year is blank')      

  return data
#______________________________________________________________________________

def getDetectionList(year,filename='default',dir='summaryStats/'):

  if filename=='default':
    filename='detectionGroups'+year+'.txt'
    
  file=open(dir+filename,'r')
    
  data=[]
  for line in file:
    fields=string.split(line,'|')
    
    if string.strip(fields[0])=='UKIRT name':
      header=[]
      for col in fields:
        header.append(string.strip(col))
    else:
      data.append({})
      for col,value in zip(header,fields):
        #print 'col',col
        #print ' value',value
        try:
          data[-1][col]=float(value)
        except:
          data[-1][col]=string.strip(value)
          
  return data
#______________________________________________________________________________

def nameSplitter(ukirtName):

#  field=ukirtName[7:9]
#  ccd=ukirtName[10:11]
# careful: two 2015 fields (5 and 6) only have 1 digit!
# more careful: 2017 field has an underscore in it!
  
  ukirtNameParts=ukirtName.split('_')
  year=ukirtName[2:6]
  if year=='2017' or year=='2018' or year=='2019':
    field=ukirtNameParts[1]+'_'+ukirtNameParts[2]
    ccd=ukirtNameParts[3]
  elif year=='2015' or year=='2016':
    field=ukirtNameParts[1]
    ccd=ukirtNameParts[2]
  else:
    print('ukirtName:',ukirtName)
    exit('ERROR: unknown year in nameSplitter')
  return field,ccd
#______________________________________________________________________________

def getNumberSources(year,photomType,field,ccd,band):

    # originally this was done by 'wc'ing in the data directories
    # now, the output of wc is saved in summaryStats/
    #  so just read that in directly

    if year=='2015' or year=='2016':
      photomfilename = photomType+'_lc_'+year[-2:]+'_'+field+'_'+ccd+'.dat'
    else:
      photomfilename = photomType+'_lc_'+year+'_'+field+'_'+ccd+'_'+band+'.dat'
        
    bruteForce=0
    if bruteForce:
      dir = options['dataDir']+'/ukirt/'+year+'/'

      if not os.access(dir+photomfilename, os.R_OK):
        print('ERROR: file missing',dir+photomfilename)
        exit('ERROR: file missing')
      else:
        os.system('wc '+dir+photomfilename+' > ztemp')
        wc=open('ztemp','r')
        os.system('rm ztemp')
        for line in wc:
            N = float(string.split(line)[0])
            #print 'N',N
        wc.close()
      return N
    
    else:
      dir = 'summaryStats/'
      countfilename = 'sourceCounts'+year+'.wcGrab'
      file = open(dir+countfilename,'r')
    
      data=[]
      for line in file:
        if photomfilename in line:
          #print 'found it!'
          cols = string.split(line)
          return float(cols[0])
      print('ERROR: didnt find the file with source counts',photomfilename)
      return 666
#______________________________________________________________________________

def readAltbandCrossCheck(year,field,ccd,
                          photomType='CASU',dir='altBandCrossCheck/'):
  import csv,string
  header=''
  data={}

  #fields,ccds = allFields(year)
  #for field in fields:
  #  for ccd in ccds:

  if 1:
      #print 'field?',field,fieldSelect
      #if field==fieldSelect  or  fieldSelect=='all':
      filename='altBandCrossCheck_'+photomType+year+'_'+field+'_'+ccd+'.txt'
      if not os.access(dir+filename, os.R_OK):
        if year=='2017' or year=='2018' or year=='2019':
          print('ERROR: altband file doesnt exist ',dir+filename)
          exit('check this!')
        return data        
      print('opening an altBand cross-check file:',filename)
      file=open(dir+filename,'r')
      file=csv.reader(file,delimiter='|')

      header=next(file)
      for i in range(len(header)):
        header[i]=header[i].strip()

      if data=={}:
        # only initialize the data dictionary for the first file
        firstheader=header

        # oops, careful: there are identical column headings in this file
        #  have to distinguish between them somehow
        # locate all duplicates and add 'alt' to their names
        #print 'verbatim header',firstheader
        cols=[]
        for i in range(len(firstheader)):
          col = firstheader[i]
          if col in cols:
            firstheader[i] = 'alt' + firstheader[i]
          else:
            cols.append(col)
        #print 'modified header',firstheader
        
        for field in header:    
          data[field]=[]
      else:
# make sure that subsequent files have the same header as previous one
#  (these two checks are redundant)
        if header!=firstheader:
          print('TROUBLE: headers dont match',header,firstheader)
          exit('ERROR: headers dont match! 1')
        for field in header:    
          if not field in data.keys():
            print('field:',field,'  keys:',data.keys())
            exit('ERROR: headers dont match! 2')

# during debugging, just read in the start of each file
      istop = 2000
      istop = 10000
      istop = -1
      icount=0
      for line in file:
        icount+=1
        if icount == istop:
          break
        
        for field,value in zip(header,line):
          if field.strip():
            #print field,'value:',value
            # float() doesn't work for the 'name' field
            try:
              data[field].append(float(value))
            except:
              data[field].append(value.strip())

  #    print '  number of data now:',len(data['name'])
  #  else:
  #    print ' skipping this field:',fieldNow

  # the final '|' creates a blank field at the end.  remove it
  if '' in data.keys():
    data.pop('')
    header.remove('')

  # convert slow lists to fast arrays
  #  be careful with the extra primary-band-only fields at the end of the file
  #  need to remove them, or arrays with be string, not float
  singleBanders = np.where(np.array(data['altN'])=='')
  Ndrop=len(singleBanders[0])
  print('  Ndata, Ndrop:',len(data['altN']),Ndrop)
  #Ndrop-=10
  for field in header:
    #print ' len ',field,len(data[field])

    if Ndrop > 0:
      data[field]=np.array(data[field])[:-Ndrop]
    #print ' len ',field,len(data[field])
    #print 'types',field,type(data[field][0])
    try:
      data[field]=np.array(data[field], dtype=np.float)
    except:
      data[field]=np.array(data[field])
    #print ' len ',field,len(data[field])
    #print 'types',field,type(data[field][0])
  #singleBanders = np.where(data['altN']=='')
  #Ndrop=len(singleBanders[0])
  #print 'Ndrop',Ndrop,len(data['altN'])

  return data
#______________________________________________________________________________

def readOverlapCrossCheck(year,photomType='CASU',dir='overlapCrossCheck/',
                          fieldSelect='all',ccdSelect='all'):
  import csv,string
  header=''
  
  fields,ccds = allFields(year)

  data={}
  for field in fields:
    if field==fieldSelect  or  fieldSelect=='all':
      filename='overlapCrossCheck_'+photomType+year+'_'+field+'.txt'
      #print 'FILENAME',filename
      if not os.access(dir+filename, os.R_OK):
        print('NOTE: file doesnt exist ',filename)
        return data        
      print(' opening an overlap cross-check file:',filename)
      file=open(dir+filename,'r')
      file=csv.reader(file,delimiter='|')

      header=next(file)
      for i in range(len(header)):
        header[i]=string.strip(header[i])

      if data=={}:
# only initialize the data dictionary for the first file
        firstheader=header

# oops, careful: there are identical column headings in this file
#  have to distinguish between them somehow
# locate all duplicates and add 'overlap' to their names
        #print 'verbatim header',firstheader
        cols=[]
        for i in range(len(firstheader)):
          col = firstheader[i]
          if col in cols:
            firstheader[i] = 'overlap' + firstheader[i]
          else:
            cols.append(col)
        #print 'modified header',firstheader
        
        for field in header:    
          data[field]=[]
      else:
# make sure that subsequent files have the same header as previous one
#  (these two checks are redundant)
        if header!=firstheader:
          print('TROUBLE: headers dont match',header,firstheader)
          exit('ERROR: headers dont match! 1')
        for field in header:    
          if not field in data.keys():
            print('field:',field,'  keys:',data.keys())
            exit('ERROR: headers dont match! 2')

# during debugging, just read in the start of each file
      istop = 2000
      istop = -1
      icount=0
      for line in file:
        icount+=1
        if icount == istop:
          break
        
# oops, actually overlap crossCheck doesnt print out the ccd info
# (but altband does now)  (modify overlap?  not a bad idea..)
        #ccdHeaderPosition=header.index('ccd')
        #ccdNow=string.strip(line[ccdHeaderPosition])
        ccdNow=ccdSelect
        #print 'ccd',ccdNow

        # special option: only select a particular chip
        if ccdNow==ccdSelect or ccdSelect=='all':
          for field,value in zip(header,line):
            if string.strip(field):
#              print field,'value:',value
# float() doesn't work for the 'name' field
              try:
                data[field].append(float(value))
              except:
                data[field].append(string.strip(value))

      #print '  number of data now:',len(data['name'])
#    else:
#      print ' skipping this field:',fieldNow

# the final '|' creates a blank field at the end.  remove it
  if '' in data.keys():
    data.pop('')
    header.remove('')
    
# convert slow lists to fast arrays
  for field in header:
    try:
      data[field]=np.array(data[field], dtype=np.float)
    except:
      data[field]=np.array(data[field])

  return data
#______________________________________________________________________________

def getCCDlocation(year,field,ccd):

  #fieldPositions.append({'year':year,'field':field,'ccd':ccd,
  #                         'RAmin':RAmin,'RAmax':RAmax,
  #                         'DECmin':DECmin,'DECmax':DECmax})

  fieldPositions = readUKIRTfields()

  for pos in fieldPositions:
      if pos['year']==year and \
         pos['field']==field and \
         pos['ccd']==ccd:

          RAmin = pos['RAmin']
          RAmax = pos['RAmax']
          DECmin = pos['DECmin']
          DECmax = pos['DECmax']
          return RAmin,RAmax,DECmin,DECmax

  exit('ERROR: year/field/ccd not found in CCD location list')          
#______________________________________________________________________________

def readUKIRTfields(filename='ukirtFieldLocations.sav',
                    dir='fieldLocations/'):

  if not os.access(dir+filename, os.R_OK):
      exit('ERROR: no file with CCD locations')

  infile = open(dir+filename,'rb')
  fieldPositions = pickle.load(infile,encoding='latin1')

  return fieldPositions
#______________________________________________________________________________


def readVVVinfo(year,options,filename='table1.tex',dir='vvv/'):

  file=open(dir+filename,'r')

  data=[]
  for line in file:
    fields=string.split(line,'&')
    
    if string.strip(fields[0])=='Tile':
      header=[]
      for col in fields:
        header.append(string.strip(col))
    else:
      data.append({})
      for col,value in zip(header,fields):
        #print 'col, value',col,value
        try:
          data[-1][col]=float(value)
        except:
          data[-1][col]=string.strip(value)
          
  return data
#______________________________________________________________________________


def readClasses(options,dir='evaluator/logData/',filename='default'):

#  import evaluator
#  import evaluator.readers 
#  from evaluator/readers import loadCurrentLog
#  execfile('evaluator/readers.py')   # hmm, doesnt import loadCurrentLog.why?
  sys.path.append('evaluator')
  from readers import loadCurrentLog

  # find the most up-to-date logfile
  options['logDir'] = dir
  if filename=='default':
    #data = loadCurrentLog(options,specifiedName=False,specifiedPerson='gb')
    data = loadCurrentLog(options,specifiedName=False)
  else:
    data = loadCurrentLog(options,specifiedName=filename)
    
  return data
#______________________________________________________________________________

def readFeaturesSimple(filename='metricList2015.txt',dir='../'):

  file=open(dir+filename,'r')

  data=[]
  for line in file:
    fields=string.split(line,'|')
    
    if string.strip(fields[0])=='name':
      header=[]
      for col in fields:
        header.append(string.strip(col))
    else:
      data.append({})
      for col,value in zip(header,fields):
        #print 'col, value',col,value
        try:
          data[-1][col]=float(value)
        except:
          data[-1][col]=string.strip(value)
          
  return data
#______________________________________________________________________________


#def readChallengeData(options, dir='dataChallenge/',
def readChallengeData(options, dir='../microlens/dataChallenge/inputData/',
                      targetSelection='All'):
  import os

  primaryband = 'W149'
  secondaryband = 'Z087'

# get the list of files in the directory
#  filelist=os.listdir(dir)
# nah. just use the list from event_info.txt

# use this factor to convert a mag err to flux err
  errorConvertFactor=0.92
# set the I band zeropoint, to convert magnitudes to Jy
  zeropointI=1.


  
  # First make a list of the event names and create an event for each
  file=open(dir+'event_info.txt','r')
  lines=file.readlines()
  file.close()

  events={} 
  for line in lines:
#  for line in lines[:10]:
    cols = line.split()
    if len(cols)!=9:
      exit('event info should have 9 columns')

    name = cols[0]
    event={'name':name,
           'RA':cols[2],
           'Dec':cols[3],
           'Distance':cols[4],
           'A_W149':cols[5],
           'sigma_A_W149':cols[6],
           'A_Z087':cols[7],
           'sigma_A_Z087':cols[8],
           'data':{}}
    events[name] = event

    
  # Now add the primary/secondary photometry to each event
  #for filename in filelist:
  for eventname in events:

    for band in [primaryband,secondaryband]:
      
      filename = 'lc/'+eventname+'_'+band+'.txt'
      print('loading a lightcurve',filename)
      file=open(dir+filename,'r')
      lines=file.readlines()
      file.close()

      times=[]
      mags=[]
      errors=[]
      for line in lines:
        cols = line.split()
        if len(cols)!=3:
          exit('should be 3 columns of data in the lightcurves')
        times.append(float(cols[0]))
        mags.append(float(cols[1]))
        errors.append(float(cols[2]))
      times=np.array(times) - options['JDadjust']
      mags=np.array(mags)
      errors=np.array(errors)
      
      photom={'time':times,
              'mag':mags,'magerr':errors}           
      photom['flux']=zeropointI * 10.**(-photom['mag']/2.5)
      photom['fluxerr']=photom['magerr']* \
                         photom['flux']*errorConvertFactor
            
      events[eventname]['data'][band]=photom

  print('# of events:',len(events.keys()))

# optionally limit the event list to some specific batch of targets
  if targetSelection!='All' and targetSelection!=[]:
    allevents=events
    print(' culling the event list',len(events))
    if len(events)<10:
      print(allevents)
    events={}
    for name in targetSelection:
#      if name:  # special check added for cases where altname is blank
# no this fails the check below.  just don't pass blank ('') in here
        #print 'targetselection name',name
        events[name]=allevents[name]
    print(' culled the event list',len(events))
    if len(events)!=len(targetSelection):
      print('NOTE: not all targets found!',len(targetSelection))
      exit('testing the target culling')

  return events
#______________________________________________________________________________

def readLog(fileName,options,useColumnThree=False):

    dir = options['logDir']
    file = open(dir+fileName,'r')

    logList = {'names':[], 'categories':[]}
    
    for line in file:
        cols = string.split(line,'|')

        LCname = string.strip(cols[0])
        if useColumnThree:
          category = string.strip(cols[2])
          # limit the final eval to just the category, not the comments
          category = category[0]
        else:
          category = string.strip(cols[1])

        if LCname in logList['names']:
            print('ERROR: a lightcurve name is in log file more than once', \
                dir+fileName,LCname)
            
        logList['names'].append(LCname)
        logList['categories'].append(category)

    file.close()

    return logList

#______________________________________________________________

def readBesties(year, minProb=0, dir='./', fileName='default'):

  if fileName=='default':
    fileName = 'besties'+year+'.txt'

  file = open(dir+fileName,'r')

  names = []
  probs = []    
  for line in file:
    cols = string.split(line)
    
    if float(cols[1]) > minProb:
      names.append(string.strip(cols[0]))
      probs.append(float(cols[1]))
      
  file.close()

  return names,probs

#______________________________________________________________

def readMLclasses(filename='classifications_finalclass_dChi300_V5_2017.txt',
                  dir='machina/'):

  file=open(dir+filename,'r')

  data=[]
  for line in file:
    fields=string.split(line,'|')
    
    if string.strip(fields[0])=='name':
      header=[]
      for col in fields:
        header.append(string.strip(col))
    else:
      data.append({})
      for col,value in zip(header,fields):
        #print 'col, value',col,value
        try:
          data[-1][col]=float(value)
        except:
          data[-1][col]=string.strip(value)
          
  return data
#______________________________________________________________________________
def loadPeriodogramResults(filename,dir='periodResults/'):
    
  file=open(dir+filename,'r')
    
  #data=[]
  data={}
  for line in file:
    fields=string.split(line,'|')
    
    if string.strip(fields[0])=='target':
      header=[]
      for col in fields:
        header.append(string.strip(col))
    else:
      #data.append({})
      thisline={}
      for col,value in zip(header,fields):
        try:
          #data[-1][col]=float(value)
          thisline[col]=float(value)
        except:
          #data[-1][col]=stringy.strip(value)
          thisline[col]=string.strip(value)
      name = thisline['target'] 
      data[name] = thisline
          
  return data
#______________________________________________________________________________

def readChallengeClasses(filename='categories.txt',
                 dir='/Users/bryden/Desktop/microlens/dataChallenge/catLogs/'):
  
  file=open(dir+filename,'r')
    
  #data={'name':[],'class':[]}
  data={}
  for line in file:
    cols = string.split(line,'|')

    #data['name'].append(string.strip(cols[0]))
    #data['class'].append(string.strip(cols[1]))
    data[string.strip(cols[0])] = string.strip(cols[1])
          
  return data

def readOverlapData(year,field,ccd,
                    saveDir='/Volumes/Data/overlapData/'):
#                    saveDir='overlapData/'):

  infileName = 'overlapData'+year+'_'+field+'_'+ccd+'.sav'
  if not os.access(saveDir+infileName, os.R_OK):
    print('infileName:',saveDir+infileName)
    exit('ERROR: overlap data is missing for this field')
  else:
    infile = open(saveDir+infileName,'r')
    dataforthisfield = {}

    Done = False
    while not Done:
      try:
        name = pickle.load(infile)
      except:
        Done = True
      if not Done:
        #print 'loading in overlap data for:',len(dataforthisfield.keys()),name
        photom = {}
        (photom['time'],photom['mag'],photom['magerr']) = pickle.load(infile)
        altphotom = {}
        (altphotom['time'],altphotom['mag'],altphotom['magerr']) = pickle.load(infile)
      
        overlapNames = pickle.load(infile)
        #print 'overlapNames',overlapNames
        overlapPhotom = []
        overlapaltPhotom = []
        if overlapNames:
          for i,ovename in enumerate(overlapNames):
            #print 'overlaploop',i
            overlapPhotom.append({})
            (overlapPhotom[i]['time'],overlapPhotom[i]['mag'],overlapPhotom[i]['magerr']) = pickle.load(infile)
            overlapaltPhotom.append({})
            (overlapaltPhotom[i]['time'],overlapaltPhotom[i]['mag'],overlapaltPhotom[i]['magerr']) = pickle.load(infile)
          
        dataforthisfield[name] = (photom,altphotom,overlapNames,overlapPhotom,overlapaltPhotom)
  
  return dataforthisfield
