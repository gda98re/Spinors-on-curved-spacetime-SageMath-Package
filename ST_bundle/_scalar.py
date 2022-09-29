from ST_bundle.tindices_object import tindices_object

from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.scalarfield import ScalarField

class _scalar(tindices_object):
    
    def __init__(self,STb,tindices_list,scalar_field = None,state = "Mutable"):
        
        if(scalar_field == None): _scalar_field = lambda : STb.M.scalar_field(0)
        else: _scalar_field = lambda : scalar_field.copy()
        super().__init__(STb,tindices_list,_scalar_field,state)
        self._chart = _scalar_field().coord_function().chart()

    @property
    def chart(self):
        return self._chart
        
    def __repr__(self):
        return f"{self.tindices} Tetrad index family of scalar fields on M with respect to the tetrad frame {self.STbundle.tframe}"
    
    def parent(self):
        return tinidces_object
    
    def __setitem__(self,pos_index,elem):
        if(not isinstance(elem,ScalarField)): raise TypeError("Can't assign something that is not a scalar field on M")
        super().__setitem__(pos_index,elem)
    
    def display(self):
        
        res = super().recursive_generator(4,super().tindices,lambda : None)
        
        lst,loop_command,_ = super().nested_loop("i",0,len(self.tindices),"","","")
        if(len(self.tindices) == 0):
            return self.comp[:]
        else:
            loop_command += f"res{lst.replace(',',']').replace('i','[i')} = self[{lst[:-1]}].expr()"
            exec(loop_command)
            return res
        
    def get_tangent_tensor(self):
        
        u = self.tindices.count("up")
        d = self.tindices.count("down")
        
        tensor_field = self.STbundle.M.tensor_field(u,d)
        tensor_field[self.STbundle.M.default_chart().frame(),[0 for k in range(0,len(self.tindices))]] = 0
        
        res = _scalar.recursive_generator(1,[""],lambda : tensor_field)
        
        super().get_tangent_tensor_alg(res,0)
        
        return res[0]
      
    def copy(self):
        res = _scalar(self.STbundle,self.tindices)
        return super().copy_alg(res)
    
    def __eq__(self,other):
        if(isinstance(other,_scalar) == False): return False
        return super().eq_alg(other)
    
    def __add__(self,other):
        res = _scalar(self.STbundle,self.tindices)
        super().sum_alg(other,res,"+")
        return res
    
    def __sub__(self,other):
        res = _scalar(self.STbundle,self.tindices)
        super().sum_alg(other,res,"-")
        return res
   
    def swap_tindices(self,pos1 = 0,pos2 = 1):
        
        res_tindices_list = self.tindices[:pos1] + [self.tindices[pos2]] + self.tindices[pos1+1:pos2] + [self.tindices[pos1]] + self.tindices[pos2+1:]
        res = _scalar(self.STbundle,res_tindices_list)
        super().swap_tindices_alg(res,pos1,pos2)
        return res
 
    def __mul__(self,other): #tensor product full (all type of indices, both tetrad and spin)
        
            if(not isinstance(other,tindices_object)): 
                res_tindices_list = self.tindices
            else:
                res_tindices_list = self.tindices + other.tindices
            
            if(isinstance(other,_scalar)):
                res = _scalar(self.STbundle,res_tindices_list)
            elif(isinstance(other,_spin_tensor)):
                res = _spin_tensor(self.STbundle,res_tindices_list,other.sindices)
            elif(isinstance(other,FreeModuleTensor)):
                res = _spin_tensor(self.STbundle,res_tindices_list,other.tensor_type())
            else:
                res = _scalar(self.STbundle,res_tindices_list)
                
            super().tensor_product_alg(other,res)
            return res
       
    def __rmul__(self,other):
        return self.__mul__(other)
    
    def __matmul__(self,other):
        return self.__mul__(other)
    
    def __rmatmul__(self,other):
        return self.__mul__(other)
        
    def ttrace(self,pos1=0,pos2=1): #trace of tindices
        
        pos1,pos2 = self.check_tindices(self,pos1,pos2,"trace")
        res_tindices_list = list(self.tindices)
        
        for k in pos1:
            res_tindices_list[k] = None
        for k in pos2:
            res_tindices_list[k] = None
        for k in range(0,2*len(pos1)):
            res_tindices_list.remove(None) 
        
        res = _scalar(self.STbundle,res_tindices_list)
        super().trace_alg(res,pos1,pos2)
        if(len(res.tindices) == 0): 
            return res.comp
        else:
            return res
          
    def tcontract(self,pos1,other,pos2):
 
        pos1,pos2 = self.check_tindices(other,pos1,pos2,"contract")
        res_tindices_list = list(self.tindices+other.tindices)
        
        for k in pos1:
            res_tindices_list[k] = None
        for k in pos2:
            res_tindices_list[k+len(self.tindices)] = None
        for k in range(0,2*len(pos1)):
            res_tindices_list.remove(None)
        
        if(not isinstance(other,tindices_object)): raise TypeError("Can't contract with objects that are not tindices_object instances") 
        if(isinstance(other,_scalar)):
            res = _scalar(self.STbundle,res_tindices_list)
        elif(isinstance(other,_spin_tensor)):
            res = _spin_tensor(self.STbundle,res_tindices_list,other.sindices)
            
        super().contract_alg(other,res,pos1,pos2)
        if(len(res.tindices) == 0): 
            return res.comp
        else:
            return res
    
    def up_down(self,pos,u_d):
        
        if(len(self.tindices) == 0 ): raise TypeError("Object passed has no tindices")
        if(self.tindices[pos] == u_d): raise TypeError(f"Index in position {pos} is already {u_d}")

        eta = self.STbundle.eta(u_d)
        res = eta*self
        res = res.swap_tindices(0,pos+2)
        res = res.ttrace(0,1)
        return res
    
    def up(self,pos=0):
        return self.up_down(pos,"up") 
    
    def down(self,pos=0):
        return self.up_down(pos,"down") 
    
    def par_der(self):
        res = self.STbundle.scalar(self.tindices+["down"])
        lst,loop_command,indent = _scalar.nested_loop("i",0,len(res.tindices),"","","")
        loop_command += f"res[{lst[:-1]}] = (self[{lst[:-4]}].derivative())(self.STbundle.tframe[{lst[-3:-1]}])"
        exec(loop_command)
        return res
    
    def tcov_der(self):
        
        Ric_coef = self.STbundle.fb_coef
        
        res = self.par_der()

        for k in range(0,len(self.tindices)):
            if(self.tindices[k] == "up"):
                par = ((self*Ric_coef).swap_tindices(k,len(self.tindices))).ttrace(len(self.tindices),len(self.tindices)+1)
                res += par
            elif(self.tindices[k] == "down"):
                par = ((self*Ric_coef).swap_tindices(k,len(self.tindices)+1)).ttrace(len(self.tindices),len(self.tindices)+1)
                res += par
        
        return res


from ST_bundle._spin_tensor import _spin_tensor

