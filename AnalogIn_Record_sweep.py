"""
   DWF Python Example
   Author:  Digilent, Inc.
   Revision: 10/17/2013

   Requires:                       
       Python 2.7, numpy, matplotlib
       python-dateutil, pyparsing
"""
import numpy as np
import matplotlib.pyplot as plt
from ctypes import *
from dwfconstants import *
import math
import time
import matplotlib.pyplot as plt
import sys

if sys.platform.startswith("win"):
    dwf = cdll.dwf
elif sys.platform.startswith("darwin"):
    dwf = cdll.LoadLibrary("/Library/Frameworks/dwf.framework/dwf")
else:
    dwf = cdll.LoadLibrary("libdwf.so")

#Frequency sweep range
startFreq = 200
endFreq = 5000
stepSize = 10
freqSpan = np.arange(startFreq,endFreq,stepSize)
amp = 2 #Set the amplitude
n = 1 #number of measurements to average over
rejectFreq = 180 #reject anything below this frequency
nSamples = 5000
sFreq = 50000

#declare ctype variables
hdwf = c_int()
sts = c_byte()
hzAcq = c_double(sFreq) #This is the sampling frequency
rgdSamples = (c_double*nSamples)()
cAvailable = c_int()
cLost = c_int()
cCorrupted = c_int()
fLost = 0
fCorrupted = 0

#print DWF version
version = create_string_buffer(16)
dwf.FDwfGetVersion(version)
print "DWF Version: "+version.value

#open device
print "Opening first device"
dwf.FDwfDeviceOpen(c_int(-1), byref(hdwf))

if hdwf.value == hdwfNone.value:
    szerr = create_string_buffer(512)
    dwf.FDwfGetLastErrorMsg(szerr)
    print szerr.value
    print "failed to open device"
    quit()

print "Preparing to read sample..."

print "Generating sine waves..."

#Initialize arrays for data storage
freqRespStorage = [[0]*len(freqSpan) for i in range(n)]
fD = [[0]*len(freqSpan) for i in range(n)]
freqDomain = [[0]*len(freqSpan) for i in range(n)]
spec = [[0]*len(freqSpan) for i in range(n)]

#set up acquisition
dwf.FDwfAnalogInChannelEnableSet(hdwf, c_int(0), c_bool(True))
dwf.FDwfAnalogInChannelRangeSet(hdwf, c_int(0), c_double(5))
dwf.FDwfAnalogInAcquisitionModeSet(hdwf, acqmodeRecord)
dwf.FDwfAnalogInFrequencySet(hdwf, hzAcq)
dwf.FDwfAnalogInRecordLengthSet(hdwf, c_double(nSamples/hzAcq.value))

for k in range(len(freqSpan)):
    if sFreq / float(nSamples) >= stepSize:
        print 'Warning!!!: Step size is less than the sampling rate, output will be undersampled'
    for l in range(n):
        dwf.FDwfAnalogOutNodeEnableSet(hdwf, c_int(0), AnalogOutNodeCarrier, c_bool(True))
        dwf.FDwfAnalogOutNodeFunctionSet(hdwf, c_int(0), AnalogOutNodeCarrier, funcSine)
        dwf.FDwfAnalogOutNodeFrequencySet(hdwf, c_int(0), AnalogOutNodeCarrier, c_double(freqSpan[k]))
        dwf.FDwfAnalogOutNodeAmplitudeSet(hdwf, c_int(0), AnalogOutNodeCarrier, c_double(amp))
        dwf.FDwfAnalogOutNodeOffsetSet(hdwf, c_int(0), AnalogOutNodeCarrier, c_double(-amp)) #negative since inverting
        dwf.FDwfAnalogOutConfigure(hdwf, c_int(0), c_bool(True))

        #wait at least 1 second(s) for the offset to stabilize
        if k==0 & l==0:
            time.sleep(1)
        else:
            time.sleep(0.05)


        #begin acquisition
        dwf.FDwfAnalogInConfigure(hdwf, c_int(0), c_int(1))
        print '   waiting to finish, frequency = '+str(freqSpan[k])+' Hz'

        cSamples = 0

        while cSamples < nSamples:
            dwf.FDwfAnalogInStatus(hdwf, c_int(1), byref(sts))
            if cSamples == 0 and (sts == DwfStateConfig or sts == DwfStatePrefill or sts == DwfStateArmed) :
                # Acquisition not yet started.
                continue

            dwf.FDwfAnalogInStatusRecord(hdwf, byref(cAvailable), byref(cLost), byref(cCorrupted))

            cSamples += cLost.value

            if cLost.value :
                fLost = 1
            if cCorrupted.value :
                fCorrupted = 1

            if cAvailable.value==0 :
                continue

            if cSamples+cAvailable.value > nSamples :
                cAvailable = c_int(nSamples-cSamples)

            # get samples
            dwf.FDwfAnalogInStatusData(hdwf, c_int(0), byref(rgdSamples, 8*cSamples), cAvailable)
            cSamples += cAvailable.value


        #print "Recording finished"
        if fLost:
            print "Samples were lost! Reduce frequency"
        if cCorrupted:
            print "Samples could be corrupted! Reduce frequency"

        f = open("record.csv", "w")
        for v in rgdSamples:
            f.write("%s\n" % v)
        f.close()

        rgpy=[0.0]*len(rgdSamples)
        for i in range(0,len(rgpy)):
            rgpy[i]=rgdSamples[i]

        freqRespStorage[l][k] = rgpy

        sample = np.asarray(freqRespStorage[l][k])
        spec[l][k] = np.fft.rfft(sample)
        dSample = hzAcq.value
        #print 'sample rate is = ' + str(dSample)
        # dSample = 1/hzAcq #Sample spacing (1/sampler rate)
        fD[l][k] = np.fft.fftfreq(len(sample), d=1/float(dSample))
        freqDomain[l][k] = np.concatenate((np.array([0]), fD[l][k][:len(sample) / 2]))

dwf.FDwfDeviceClose(hdwf)

#plt.figure(1)
#plt.plot(rgpy)
#plt.show()

#Sum signals at same frequency
m = np.asarray(spec)
q = sum(m)

#Compute the spectrum
maximum = [[0]*len(freqSpan) for x in range(n)]
for k in range(len(freqSpan)):
    #normSamples = np.asarray([float(abs(x)) / n for x in q[k]])  # normalize summed signals
    for l in range(n):
        normSamples = np.asarray([abs(x) for x in m[l][k]])  # normalize summed signals
        rejectIndex = int(rejectFreq/(sFreq/float(nSamples)))
        lowerIndex = int(max([rejectIndex,(freqSpan[k]/(sFreq/float(nSamples)))-2]))
        upperIndex = int(freqSpan[k]/(sFreq/float(nSamples)))+10
        maximum[l][k] = max(normSamples[lowerIndex:upperIndex])
        #summation = m_r + m_i
        #absolute = abs(summation)
        print 'Max absolute at '+str(freqSpan[k])+' Hz sample '+str(l)+' = ' + str(maximum[l][k])

        plt.figure(1)
        plt.plot(freqDomain[l][k],normSamples,(0.5, 0.5, 0.5))
        #plt.plot(freqDomain,abs(spec.real),'-r',freqDomain,abs(spec.imag),'-k')
        plt.show()
        #plt.show(block=False)
        #plt.pause(5)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amp (???)')
        plt.title('FFT of received signal')
        axes = plt.gca()
        axes.set_xlim([0,5000])
        #axes.set_ylim([0,1000])

SPL_Out = [float(x)/n for x in sum(np.asarray(maximum))]
plt.figure(2)
plt.plot(freqSpan,SPL_Out)