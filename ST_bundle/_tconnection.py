from ST_bundle._spin_tensor import _spin_tensor

from sage.manifolds.differentiable.bundle_connection import BundleConnection
from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.differentiable.vectorfield import VectorField

class _tconnection(BundleConnection):

    def __init__(self,STb):
        
        BundleConnection.__init__(self,STb.fbundle,"\\nabla_t")
        
        self._STb = STb
        
        for i in range(0,4):
            for j in range(0,4):
                self.set_connection_form(i,j)[:] = STb.connection.connection_form(j,i,STb.tframe)
     
    @property   
    def STbundle(self):
        return self._STb
    
    @property
    def connection_coef(self):
        return self.STbundle.fb_coef
    
    def __call__(self,tangent_section,section):
        
        #D(V,F) #F section of frame bundle
        if(isinstance(tangent_section,VectorField) and section.parent() == self.STbundle.fbundle.section_module()):    
            return BundleConnection.__call__(self,tangent_section,section)
        
        #D(e_{(a)},F)
        if(tangent_section == self.STbundle.tframe and section.parent() == self.STbundle.fbundle.section_module()):
            res = self.STbundle.scalar(["down"])
            for i in range(0,4):
                res[i] = BundleConnection.__call__(self,tangent_section[i],section)
            return res
        
        #convert a tensor on the spin bundle to a tindices object with no tindices or a tensor on the frame bundle to a scalar with tindces objects                                                        
        if(isinstance(section,FreeModuleTensor) and section.base_module() == self.STbundle.fbundle.section_module()): 
            tindices_list = ["up" for k in range(0,section.tensor_type()[0])]+["down" for k in range(0,section.tensor_type()[1])]
            section = self.STbundle.set_tscalar_from_Components(tindices_list,section.comp())
        elif(isinstance(section,FreeModuleTensor) and section.base_module() == self.STbundle.sbundle.section_module()):
            section = self.STbundle.spin_tensor([],section.tensor_type())
                                                                 
        res = section.tcov_der()                                                         
        #D(V,T_{(a)...}^{(b)...})                                                         
        if(isinstance(tangent_section,VectorField)):
            comp = self.STbundle.set_tscalar_from_Components(["up"],tangent_section.comp(self.STbundle.tframe))
            res = comp.tcontract(0,res,len(res.tindices)-1)
            return res
                
        #D(e_{(c)},T_{(a)...}^{(b)...})  note: index (c) is the last in the object contructed                                                      
        elif(tangent_section == self.STbundle.tframe):
            return res
