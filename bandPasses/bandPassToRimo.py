__DESCRIPTION__="""Convert bandpasses per diode into a fits file formatted as the RIMO"""
 
import numpy as np
import sys
import os
import pyfits
from SmartTable import csv_table

def diode_arms_association(a) :
   if type(a) == type('') :
      code = a[0]
   elif type(a) == type(0) :
      code = 0
   else :
      return
   if a==0 or a=='0' or a=='M' or a=='m' or a=='Y' or a=='y' : return {'diodes':['r00','r01'],'arm':'M','alt_arm':'Y'}
   if a==1 or a=='1' or a=='S' or a=='s' or a=='X' or a=='x' : return {'diodes':['r10','r11'],'arm':'S','alt_arm':'X'}
   return

class _baseBP :
   def __init__(self) : pass
   def get_dat(self) : pass
   def __len__(self) : 
      try :
         return len(self.Freq)
      except :
         return 0
   def pickle(self,picklefile) :
      import pickle
      try :
         pickle.dump(self.__dict__,open(picklefile,'w'))
      except :
         print "Impossible to create ",picklefile
   def load(self,picklefile) :
      import pickle
      try :
         self.__dict__=pickle.load(open(picklefile,'r'))
      except :
         print "Impossible to open ",picklefile
   def toFits(self) : pass
   def diodes(self) :
      """returns list of diodes"""
      return ['r00','r01','r10','r11']
   def channel(self) :
      if self.fh==None :
         return 0
      if self.fh < 18 :
         return 0
      elif self.fh < 24 : 
         return 70 
      elif self.fh <27 : 
         return 44
      elif self.fh <29 : 
         return 30
      else : 
         return 0
      
class qucsBP(_baseBP) :
   def __init__(self,fh,diode,path='LFI_bandpasses_QUCS_20091109',template='LFI_RCA%d_%dGHz_%s-%s.dat',filename=None):
      """Object to handle dat files with qucs bp"""
      self.doc=['bandpasses QUCS model from 2009 11 09 issue']
      self.fh=fh
      self.diode=diode
      self.filename = None
      self.Freq=None
      if fh==None or diode==None : return
      if filename != None :
         self.filename=filename
      else :
         das=diode_arms_association(diode[1])
         self.filename=path+'/'+template%(int(fh),int(self.channel()),das['arm'],diode[1:3])
      self.get_dat() 
   def get_dat(self,fsep='[\s\t,]+',useRegEx=True) :
      import re
      try :
         f=open(self.filename,'r') 
      except :
         print "impossible to open ",self.filenae
         return
      self.Freq=[]
      self.Trans=[]
      for l in f :
         if l.strip() != "" :
            if useRegEx :
               a=np.array(re.split(fsep,l.strip()),dtype='float')
            else :
               a=np.array(l.strip().split(fsep),dtype='float')
            self.Freq.append(a[0])
            self.Trans.append(a[1])
      for l in ['Freq','Trans'] : self.__dict__[l]=np.array(self.__dict__[l])
      self.step=self.Freq[1]-self.Freq[0]
      self.Trans *= 1/(self.step*self.Trans.sum())
      f.close()
   def toFits(self) :
      import pyfits
      from collections import OrderedDict
      import time
      import dict2fits
      hdulist=[]
      das = diode_arms_association(self.diode[1])
      arm=diode_arms_association(self.diode[1])['arm']
      alt_arm=diode_arms_association(self.diode[1])['alt_arm']
      print self.diode,' ',self.diode[1],' ',arm,' ',alt_arm,self.step,
      extname='BANDPASS_%03d-%2d%s-%s-%s'%(self.channel(),int(self.fh),arm,alt_arm,self.diode)
      print extname
      tb=OrderedDict()
      tb['WAVENUMBER']=self.Freq
      tb['TRANSMISSION']=self.Trans
      tb['UNCERTAINTY']=np.zeros(len(self))
      tb['FLAG']=np.repeat('F',len(self))
      hdulist=dict2fits.Table(tb,table_name=extname)
      hdulist.header.update('tunit1','GHz')
      hdulist.header.update('tunit2','n/a')
      hdulist.header.update('tunit3','n/a')
      hdulist.header.update('tunit4','n/a')
      hdulist.header.update('FILE',self.filename)
      hdulist.header.update('fh',self.fh,'the feed horn')
      hdulist.header.update('arm',arm,'the polarization arm')
      hdulist.header.update('altarm',alt_arm,'alternate name for arm')
      hdulist.header.update('diode',self.diode,'the diode')
      hdulist.header.update('step',self.step,'[GHz] freq step')
      hdulist.header.update('created',time.asctime(),'creation date')
      hdulist.header.add_comment('diodes 00/01 - M - X')
      hdulist.header.add_comment('diodes 10/11 - S - Y')
      hdulist.add_checksum()
      return hdulist

class diodeWeights(csv_table) :
   def __init__(self,*arg) :
      csv_table.__init__(self,None)
      csv_table.__keys__={}
      if len(arg)== 0 : return
      self.get_ascii(arg[0])
   def get_ascii(self,filename,fsep='[\s\t,]+',useRegEx=True) :
      import re
      import numpy as np
      try :
         f=open(filename,'r') 
      except :
         print "impossible to open ",self.filename
         return
      self.Name=[]
      self.w=[]
      self.fh=[]
      self.polarization=[]
      self.diode=[]
      for l in f :
         if l.strip() != "" :
            if useRegEx :
               a=re.split(fsep,l.strip())
            else :
               a=l.strip().split(fsep)
            self.Name.append(a[0])
            self.w.append(float(a[1]))
            self.fh.append(int(a[0][3:5]))
            self.polarization.append(a[0][5])
            self.diode.append(a[0][7:])
      f.close()
      for k in self.__dict__.keys() : self.__dict__[k]=np.array(self.__dict__[k])
      self.filename=filename
   def __len__(self) : return len(self.Name)

class measuredBP(_baseBP) :
   def __init__(self,fh,path='measures',template='LFI%d_%dGHz_band_meas.dat',filename=None):
      """Object to handle dat files with measured bp"""
      self.doc=['bandpasses from measures in 2007 data']
      self.fh=fh
      self.filename = None
      self.Freq=None
      if filename != None :
         self.filename=filename
         #self.channel=0
      else :
         if fh < 18 :
            return 
         elif fh < 24 : 
            freq=70 
         elif fh <27 : 
            freq=44
         elif fh <29 : 
            freq=30
         else : 
            return 
         #self.channel=int(freq)
         self.filename=path+'/'+template%(int(fh),int(self.channel()))
      self.get_dat() 
   def get_dat(self,fsep='[\s\t,]+',useRegEx=True) :
      import re
      try :
         f=open(self.filename,'r') 
      except :
         print "impossible to open ",self.filenae
         return
      self.Freq=[]
      self.r00=[]
      self.r01=[]
      self.r10=[]
      self.r11=[]
      for l in f :
         if l.strip() != "" :
            if useRegEx :
               a=np.array(re.split(fsep,l.strip()),dtype='float')
            else :
               a=np.array(l.strip().split(fsep),dtype='float')
            #a=np.array(l.strip().split('\t'),dtype='float')
            self.Freq.append(a[0])
            self.r00.append(a[1])
            self.r01.append(a[2])
            self.r10.append(a[3])
            self.r11.append(a[4])
      for l in ['Freq','r00','r01','r10','r11'] : self.__dict__[l]=np.array(self.__dict__[l])
      self.step=self.Freq[1]-self.Freq[0]
      for l in ['r00','r01','r10','r11'] : self.__dict__[l] *= 1/(self.step*self.__dict__[l].sum())
      f.close()
   def toFits(self) :
      import pyfits
      from collections import OrderedDict
      hdulist=[]
      for diode in self.diodes() :
         das = diode_arms_association(diode[1])
         arm=diode_arms_association(diode[1])['arm']
         alt_arm=diode_arms_association(diode[1])['alt_arm']
         print diode,' ',diode[1],' ',arm,' ',alt_arm,self.step,
         extname='BANDPASS_%03d-%2d%s-%s-%s'%(self.channel(),self.fh,arm,alt_arm,diode)
         print extname
         tb=OrderedDict()
         tb['WAVENUMBER']=self.Freq
         tb['TRANSMISSION']=self.__dict__[diode]
         tb['UNCERTAINTY']=np.zeros(len(self))
         tb['FLAG']=np.repeat('F',len(self))
         hdulist.append(dict2fits.Table(tb,table_name=extname))
         hdulist[-1].header.update('tunit1','GHz')
         hdulist[-1].header.update('tunit2','n/a')
         hdulist[-1].header.update('tunit3','n/a')
         hdulist[-1].header.update('tunit4','n/a')
         hdulist[-1].header.update('FILE',self.filename)
         hdulist[-1].header.update('COLUMN',diode)
         hdulist[-1].header.update('fh',fh,'the feed horn')
         hdulist[-1].header.update('arm',arm,'the polarization arm')
         hdulist[-1].header.update('altarm',alt_arm,'alternate name for arm')
         hdulist[-1].header.update('diode',diode,'the diode')
         hdulist[-1].header.update('step',self.step,'[GHz] freq step')
         hdulist[-1].header.update('created',time.asctime(),'creation date')
         hdulist[-1].header.add_comment('diodes 00/01 - M - X')
         hdulist[-1].header.add_comment('diodes 10/11 - S - Y')
      for k in range(len(hdulist)) : 
         try :
            hdulist[k].add_checksum()
            print "Checksum for hdu ",k
         except :
            pass
      return hdulist

if __name__=='__main__' :
   import time
   import pyfits
   import dict2fits
   from dict2fits import TableInfos
   from collections import OrderedDict
   DiodeWeights=diodeWeights('LFI_bandpasses_QUCS_20091109/weights_march15_2010.txt')
   
   diode_radiometer={'0':'MY','1':'SX'}
   arm_alias={'x':'S','y':'M'}
   alias_arm={'S':'x','M':'y'}
   
   codeName='bandPassToRimo.py'
   codeVersion='0.0'
   codeDate='2013 Jan 24'
   codeAuthor='M.Maris'
   
   # measured bandpasses
   FITSFILE ='LFI_bandpasses_measured_by_diode.fits'
   hdulist=[]
   # creates primary header
   hdulist.append(pyfits.PrimaryHDU())
   hdulist[-1].header.update('CName',codeName)
   hdulist[-1].header.update('CVersion',codeVersion)
   hdulist[-1].header.update('CAuthor',codeAuthor)
   hdulist[-1].header.update('CDate',codeDate)
   hdulist[-1].header.update('Created',time.asctime())
   for k in measuredBP(None).doc :   hdulist[-1].header.add_comment(k)
   hdulist[-1].header.update('extname','bp_meas_fits')
   hdulist[-1].add_checksum()
   print "Checksum for hdu ",0
   # creates tables for each radiometer/detector
   for fh in [18,19,20,21,22,23,24,25,26,27,28] :
      for k in measuredBP(fh).toFits() : 
         hdulist.append(k)
   thdulist=pyfits.HDUList(hdulist)
   thdulist.writeto(FITSFILE,clobber=True)
   
   # qucs bandpasses
   FITSFILE ='LFI_bandpasses_qucs_by_diode.fits'
   hdulist=[]
   # creates primary header
   hdulist.append(pyfits.PrimaryHDU())
   hdulist[-1].header.update('CName',codeName)
   hdulist[-1].header.update('CVersion',codeVersion)
   hdulist[-1].header.update('CAuthor',codeAuthor)
   hdulist[-1].header.update('CDate',codeDate)
   hdulist[-1].header.update('CREATED',time.asctime())
   for k in qucsBP(None,None).doc :   hdulist[-1].header.add_comment(k)
   hdulist[-1].header.update('extname','bp_qucs_fits')
   hdulist[-1].add_checksum()
   print "Checksum for hdu ",0
   # creates tables for each radiometer/detector
   for fh in [18,19,20,21,22,23,24,25,26,27,28] :
      for diode in qucsBP(None,None).diodes():
         print fh,diode
         hdulist.append(qucsBP(fh,diode).toFits())
   thdulist=pyfits.HDUList(hdulist)
   thdulist.writeto(FITSFILE,clobber=True)
   
