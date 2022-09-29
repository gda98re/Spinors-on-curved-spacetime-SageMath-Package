from ST_bundle._spin_tensor import _spin_tensor
from ST_bundle._scalar import _scalar

from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.differentiable.vectorfield import VectorField
from sage.manifolds.differentiable.bundle_connection import BundleConnection

class _totconnection:

    def __init__(self,STb,A = None):
                
        self._sconnection = STb.sconnection(A)
        self._tconnection = STb.tconnection()

    def __repr__(self):
        return f"Total connection"      

    def display(self,typ):

        if(typ == "t"):
            return self.tconnection.display()
        elif(typ == "s"):
            return self.sconnection.display()
        else:
            raise TypeError("Argument must be either 't' or 's'")


    @property
    def STbundle(self):
        return self.sconnection.STbundle

    @property    
    def sconnection(self):
        return self._sconnection

    @property    
    def tconnection(self):
        return self._tconnection

    @property    
    def sconnection_coef(self):
        return self.sconnection.connection_coef

    @property    
    def tconnection_coef(self):
        return self.tconnection.connection_coef
    
    def __call__(self,tangent_section,section):
        
        if(isinstance(section,_spin_tensor)):
            if(len(section.tindices) == 0):
                return (self.sconnection)(tangent_section,section) #oggetto con soli sindices -> sconnection
            else:
                return section.totcov_der(self.sconnection_coef) #oggetto con tinidces e sindices -> total connection
        elif(isinstance(section,_scalar)):
            return (self.tconnection)(tangent_section,section) #oggetto con soli tinidces -> tconnection
            
        if(section.base_module() == self.STbundle.fbundle.section_module()): #tensore di sage del frame bundle -> tconnection
            return (self.tconnection)(tangent_section,section)
            
        if(section.base_module() == self.STbundle.sbundle.section_module()): #tensore di sage dello spinor bundle -> sconnection
            return (self.sconnection)(tangent_section,section)
            
        
