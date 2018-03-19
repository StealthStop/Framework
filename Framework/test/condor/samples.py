from ctypes import cdll
from ctypes import c_char_p
from ctypes import c_double
from ctypes import POINTER 

class SampleCollection:
    def __init__(self):
        self.lib = cdll.LoadLibrary('../obj/samplesModule.so')
        self.obj = self.lib.SC_new()
        self.lib.SC_samples.restype = POINTER(c_char_p)
        self.lib.SC_samples_names.restype = POINTER(c_char_p)
        self.lib.SC_samplecollection_names.restype = POINTER(c_char_p)
        self.lib.SC_samplecollection_lumis.restype = POINTER(c_double)
        self.lib.SC_fixed_lumi.restype = c_double

    def nSamples(self, name):
        return self.lib.SC_samples_size(self.obj, name)

    def getFixedLumi(self):
        return self.lib.SC_fixed_lumi()

    def sampleList(self, name):
        names = self.lib.SC_samples(self.obj, name)
        files = self.lib.SC_samples_names(self.obj, name)
        list = [(names[i],files[i]) for i in xrange(self.lib.SC_samples_size(self.obj, name))]
        return list

    def sampleCollectionList(self):
        names = self.lib.SC_samplecollection_names(self.obj)
        list = [names[i] for i in xrange(self.lib.SC_samplecollection_size(self.obj))]
        return list

    def sampleCollectionLumiList(self):
        names = self.lib.SC_samplecollection_names(self.obj)
        lumis = self.lib.SC_samplecollection_lumis(self.obj)
        list = [(names[i],lumis[i]) for i in xrange(self.lib.SC_samplecollection_size(self.obj))]
        return dict(list)

