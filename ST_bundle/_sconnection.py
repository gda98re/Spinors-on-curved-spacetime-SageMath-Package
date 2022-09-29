from ST_bundle._spin_tensor import _spin_tensor

from sage.manifolds.differentiable.bundle_connection import BundleConnection
from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.differentiable.vectorfield import VectorField

class _sconnection(BundleConnection):

    def __init__(self,STb,A = None):
        
        BundleConnection.__init__(self,STb.sbundle,"\\nabla_s")
        
        self._STb = STb
        
        if(A == None): #setting A to zero
            A = 0*STb.tcoframe[0]
        
        self._A_one_form = A
        
        self._A_tetrad = STb.set_tscalar_from_Components(["down"],A.comp(STb.tframe))
        
        Id = STb.Identity([])
        
        self._connection_coef = STb.spinor_coef+Id*self.A_tetrad
        
        connection_one_forms = STb.spinor_coef.get_tangent_tensor()
        
        
        if(A == None):
            for i in range(0,4):
                for j in range(0,4):
                    super(_sconnection,self).set_connection_form(i,j)[:] = connection_one_forms[j][i]
        else:
            for i in range(0,4):
                for j in range(0,4):
                    super(_sconnection,self).set_connection_form(i,j)[:] = connection_one_forms[j][i] + Id.comp[j,i]*A
        
    @property
    def STbundle(self):
        return self._STb
    
    @property
    def A_one_form(self):
        return self._A_one_form
  
    @property
    def A_tetrad(self):
        return self._A_tetrad
    
    @property
    def connection_coef(self):
        return self._connection_coef
    
    def __call__(self,tangent_section,s_section):
        
        #D(V,psi)
        if(isinstance(tangent_section,VectorField) and s_section.parent() == self.STbundle.sbundle.section_module()): 
            return BundleConnection.__call__(self,tangent_section,s_section)
        
        #D(e_{(a)},psi)
        if(tangent_section == self.STbundle.tframe and s_section.parent() == self.STbundle.sbundle.section_module()):
            res = self.STbundle.spin_tensor(["down"],(1,0))
            for i in range(0,4):
                res[i] = BundleConnection.__call__(self,tangent_section[i],s_section)
            return res
        
        #convert a tensor on the spin bundle to a tindices object with no tindices                                                         
        if(isinstance(s_section,FreeModuleTensor) and s_section.base_module() == self.STbundle.sbundle.section_module()):
            s_section = self.STbundle.spin_tensor([],s_section.tensor_type())
                                                                 
        res = s_section.scov_der(self.connection_coef)                                                         
        #D(V,T_{(a)...}^{(b)...})                                                         
        if(isinstance(tangent_section,VectorField)):
            comp = self.STbundle.set_tscalar_from_Components(["up"],tangent_section.comp(self.STbundle.tframe))
            res = comp.tcontract(0,res,len(res.tindices)-1)
            return res
                
        #D(e_{(c)},T_{(a)...}^{(b)...})  note: index (c) is the last in the object contructed                                                      
        elif(tangent_section == self.STbundle.tframe):
            return res
            
        
