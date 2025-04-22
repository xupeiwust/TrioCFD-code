
import numpy as np
from cmath import sin
from math import *
import pickle


class Field(object):
    """ Field object is used for scalar, vectorial or tensorial quantity in order to operate between them all
     the operation we want, for instance product of a vector by a matrix or tensorial product between a vector and
     a 3D tensor ans so on."""

    def __init__(self, ax=np.array([]), npa=np.array([])):
        """ 
        initialization of a field object
        @param a np.array([]) for the x-axis : ax
        @param a np.array([]) for the y-axis : npa    
        """
        self.ax = ax
        self.npa = npa       
   
    
    def setax(self, ax):
        self.ax=ax
    
    def setnpa(self, npa):
        self.npa=npa
                
    def __add__(self, other):
        """ 
        @param re-define of the addition operator +
        @self a field
        @other a field
        """ 
        if isinstance(other,Field):
            return Field(ax=self.ax, npa=self.npa+other.npa)
        else:
            return Field(ax=self.ax, npa=self.npa+other)

    def __radd__(self, other):
        return Field(ax=self.ax, npa=self.npa + other)
    
    
    def __sub__(self, other):
        """ 
        @param re-define of the addition operator -
        """
        if isinstance(other,Field):
            return Field(ax=self.ax, npa=self.npa-other.npa)
        else:
            return Field(ax=self.ax, npa=self.npa-other)  
        
        
    def __rsub__(self, other):
        """ 
        @param re-define of the product operator *
        """
        return Field(ax=self.ax, npa=other-self.npa)  
        

    def __mul__(self, other):
        """ 
        @param re-define of the product operator *
        """
        if isinstance(other,Field):
            return Field(ax=self.ax, npa=self.npa*other.npa)
        else:
            return Field(ax=self.ax, npa=self.npa*other) 
        
    def __rmul__(self, other):
        """ 
        @param re-define of the product operator *
        """
        return Field(ax=self.ax, npa=self.npa*other) 
    
    
    def __truediv__(self, other):
        """ 
        @param re-define of the product operator *
        """
        if isinstance(other,Field):
            return Field(ax=self.ax, npa=self.npa/other.npa)
        else:
            return Field(ax=self.ax, npa=self.npa/other)  
        
        
    def __rtruediv__(self, other):
        """ 
        @param re-define of the product operator *
        """
        return Field(ax=self.ax, npa=other/self.npa)
        
        
    def __pow__(self, other):
        """ 
        @param re-define of the product operator *
        """
        return Field(ax=self.ax, npa=self.npa**other)         
        
        
    def __str__(self):
        return "\n PRINT FIELD \n Y_AXIS : {0} \n X_AXIS : {1}".format(self.npa,self.ax)   
    
    
    def __len__(self):
        return len(self.npa)
    
    def __getitem__(self, key):
        return self.npa[key]

    def __setitem__(self, key, value):
        self.npa[key] = value
    
    
    def mean(self):
        """ 
        Return the mean value of self
        @param a Field (scalar Field)
        """         
        return np.mean(self.npa) 
        
    def max(self):   
        return np.max(self.npa) 
        
    def min(self):   
        return np.min(self.npa) 
    
    def sqrt(self):
        """ 
        Return the square root of self
        @param a Field (scalar Field)
        """ 
        return Field(ax=self.ax, npa=np.sqrt(self.npa))  
    
    def exp(self):
        return Field(ax=self.ax, npa=np.exp(self.npa))  

    
    def abs(self):
        return Field(ax=self.ax, npa=np.abs(self.npa))    
    
    
    def gradient(self):
        return Field(ax=self.ax, npa=np.gradient(self.npa, self.ax))

    def derivatePlus(self):
        return Field(ax=self.ax[:-1], npa=(self.npa[1:]-self.npa[:-1])/(self.ax[1:]-self.ax[:-1]))

    def derivateMoins(self):
        return Field(ax=self.ax[1:], npa=(self.npa[1:]-self.npa[:-1])/(self.ax[1:]-self.ax[:-1]))

    def SortFromXvalue(self):
        ax_sorted=sorted(self.ax)
        npa_sorted=[self.npa for self.ax, self.npa in sorted(zip(self.ax,self.npa))]
        return Field(ax=ax_sorted, npa=npa_sorted)


    def CloserValueForXequal(self, target):
        return self.npa[min(range(len(self.ax)), key=lambda i: abs(self.ax[i]-target))]

    def CloserIndexForYequal(self, target):
        return min(range(len(self.ax)), key=lambda i: abs(self.npa[i]-target))

    def CloserIndexForXequal(self, target):
        return min(range(len(self.ax)), key=lambda i: abs(self.ax[i]-target))

    def integrale(self,amin=-1E7,amax=1E7, normalisation=True):
        integral = 0.
        norme = 0.
        for i in range(len(self.ax)-1) :
            if (self.ax[i]>amin and self.ax[i]<amax):
                integral = integral + self.npa[i]*(self.ax[i+1]-self.ax[i])
                norme = norme + self.ax[i+1]-self.ax[i]
        if normalisation : 
            return integral/norme
        else:
            return integral
        
    def ValeurCumulee(self,amin=-1E7,amax=1E7, normalisation=True):
        
        results=self*0.
        for i in range(len(self.ax)):
                results.npa[i]=self.integrale(amin,self.ax[i],normalisation)
    
        return results       
        
    
    def integrale_glissante(self,amin=-1E7, normalisation=True):
        integral = self * 0.
        norme = 0.
        for j in range(len(integral.ax)):
            integral.npa[j] = 0.
            norme=0.
            for i in range(len(self.ax)-1) :
                if (self.ax[i]>amin and self.ax[i]<integral.ax[j]):
                    integral.npa[j] = integral.npa[j] + self.npa[i]*(self.ax[i+1]-self.ax[i])
                    norme = norme + self.ax[i+1]-self.ax[i]
            if normalisation : 
                integral.npa[j] = integral.npa[j]/norme
 
        return integral
    
    def ProjectionOnAnOtherField(self, field_ref):
        resu=Field(ax=field_ref.ax*1., npa=field_ref.npa*1.)
        for i in range(len(field_ref.ax)):
            resu[i]=self.CloserValueForXequal(field_ref.ax[i])
        return resu
    
    
    def append(self, other):
        return Field(ax=np.append(self.ax, other.ax), npa=np.append(self.npa, other.npa))    
        
           
    def SupprimerSautPerio(self, amplitude_max):
        decallage_cumule = self * 0.
        for i in range(len(self.npa)):
    	    if (i==0 or i>=len(self.npa)-1):
    	        pass
    	    elif (abs(self.npa[i]-self.npa[i-1])>amplitude_max):
    	    	decallage_cumule[i:] -= self.npa[i]-self.npa[i-1]-self.npa[i+1]+self.npa[i]
        for i in range(len(self.ax)):
            self.npa[i]+=decallage_cumule[i] 
        return 


    @classmethod
    def LoadFromFile(cls, filename, Xcol, Ycol, separator=None):
        f = open(filename,"r")
        datalines = f.readlines()
        f.close()
        Y=[]
        X=[]
        for line in datalines:
            if (line[0] != "@") and (line[0] != "#") and ( len(line.split(separator)) != 0 ):
                datalist = line.split(separator)
                try:
                    try:
                        X.append(float(datalist[Xcol]))
                    except:
                        X.append(float("nan"))
                except:
                    try:
                        X.append(datalist[Xcol])  
                    except:
                        X.append(float("nan"))
                try:
                    try:
                        Y.append(float(datalist[Ycol]))
                    except:
                        Y.append(float("nan"))                       
                except:
                    try:
                        Y.append(datalist[Ycol])   
                    except:
                        Y.append(float("nan")) 
        return Field(ax=np.asarray(X),npa=np.asarray(Y))
    
    @classmethod
    def LoadFromNCFDprofile(cls, filename, Xcol):
        f = open(filename,"r")
        datalines = f.readlines()
        f.close()
        dictionnaire={}
        for line in datalines:
            var=line.split()
            if (line.split()[0] == "#COLUMN_TITLES:"):
                var.remove("#COLUMN_TITLES:")
                while('|' in var):
                    var.remove('|')
                for i in range(len(var)) :
                    dictionnaire[var[i]]=cls.LoadFromFile(filename, Xcol, i)
        return dictionnaire
        
    @classmethod
    def LoadFromTrioIJKStats(cls, filename, only_n_first_column=1e5):
        f = open(filename,"r")
        datalines = f.readlines()
        f.close()
        dictionnaire={}
        for line in datalines:
            var=line.split()
            if (line.split()[0] == "#"):
                var.remove("#")
                for i in range(len(var)) :
                    if(i<only_n_first_column):
                    	dictionnaire[var[i]]=cls.LoadFromFile(filename, 0, i)
        return dictionnaire
        
    @classmethod
    def LoadFromTrioIJKStats_brute(cls, filename, only_n_first_column=1e5):
        f = open(filename,"r")
        datalines = f.readlines()
        f.close()
        dictionnaire={}
        for line in datalines:
            var=line.split()
            if (line.split()[0] == "#" and line.split()[1]=="colonne" ):
                colonne=int(line.split()[2])-1
                var=line.split()[4]
                if (colonne<only_n_first_column):
                    dictionnaire[var]=cls.LoadFromFile(filename, 0, colonne)
        return dictionnaire


    @classmethod
    def LoadScalarFromFileWithPattern(cls, filename, Xcol, pattern):
        f = open(filename,"r")
        datalines = f.readlines()
        f.close()
        for line in datalines:
            if (line[0] == pattern) :
                datalist = line.split()
        return float(datalist[Xcol])             


    @classmethod    
    def LoadFromEvol(cls, filename):
        Y=[]
        X=[]
        try:
            f = open(filename,"r")
        except:
            print (filename, "does not exist")
            return Field(ax=np.asarray(X),npa=np.asarray(Y)) 
        datalines = f.readlines()
        f.close()

        for line in datalines:
            if (line[0] != "@") and (line[0] != "#") and ( len(line.split()) != 0 ):
                datalist = line.split()
                try:
                    X.append(float(datalist[0]))
                except:
                    X.append(datalist[0])                    
                try:
                    Y.append(float(datalist[1]))
                except:
                    Y.append(datalist[1])                    
        return Field(ax=np.asarray(X),npa=np.asarray(Y))  
        




