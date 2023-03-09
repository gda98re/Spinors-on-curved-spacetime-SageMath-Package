from ST_bundle.tindices_object import tindices_object

from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.scalarfield import ScalarField

class _scalar(tindices_object):
    
    def __init__(self,STb,tindices_list,scalar_field = None,state = "Mutable"):
        super().__init__(STb,tindices_list,ScalarField,state)
        if(scalar_field != None):
            _scalar_field = lambda : scalar_field.copy()
            lis = _scalar.recursive_generator(4,tindices_list,_scalar_field)
            if(self._nid == 0):
                self.comp = lis
            else:
                self[:] = lis
        
        
    def __repr__(self):
        return f"{self.tindices} Tetrad index family of scalar fields on M with respect to the tetrad frame {self.STbundle.tframe}"
    
    def parent(self):
        return tindices_object
                  
    def copy(self):
        res = _scalar(self.STbundle,self.tindices)
        return super()._copy_(res)
    
    def swap_adjacent_tindices(self, pos1, pos2, pos3,ref = "copy"):
        
        res_tindices_list = self.tindices[:pos1] + self.tindices[pos2:pos3] + self.tindices[pos1:pos2] + self.tindices[pos3:]
        res = _scalar(self.STbundle,res_tindices_list)
        return super()._swap_adjacent_tindices_(res,pos1, pos2, pos3,ref)

    def swap_tindices(self, pos1 = 0, pos2 = 1,ref = "copy"):
        
        res_tindices_list = self.tindices[:pos1] + [self.tindices[pos2]] + self.tindices[pos1+1:pos2] + [self.tindices[pos1]] + self.tindices[pos2+1:]
        res = _scalar(self.STbundle,res_tindices_list)
        return super()._swap_tindices_(res,pos1,pos2,ref)

    
    def __pos__(self): #+
        return self.copy()
    
    def __neg__(self): #-
        res = self.copy()
        return super()._neg_(res)

    
    def __add__(self,other):
        res = self.copy()
        return super()._add_(other,res)

    

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
                
            return super()._tensor_product_(other,res) # (*)

    def __rmul__(self,other):
        return self.__mul__(other)   
    
    def __matmul__(self,other):
        return self.__mul__(other)

    def __rmatmul__(self,other):
        return self.__mul__(other)
    
    def __truediv__(self,other):
        res = self.copy()
        return super()._truediv_(other,res)
    
        
    def trace(self,*args): #trace of tindices
        
        nargs = len(args)
        if(nargs == 0): 
            if(self.tindices == ["up","down"] or self.tindices == ["down","up"]):
                tpos1 = 0
                tpos2 = 1
            else:
                raise TypeError("At least two indices must be provided")
        elif(nargs == 1):
            raise TypeError("At least two indices must be provided")
        elif(nargs == 2):
            tpos1 = args[0]
            tpos2 = args[1]
        else:
            raise TypeError("Too many arguments") 
        
        tpos1,tpos2 = self.check_tindices(self,tpos1,tpos2,"trace")
        res_tindices_list = list(self.tindices)
        
        for k in tpos1:
            res_tindices_list[k] = None
        for k in tpos2:
            res_tindices_list[k] = None
        for k in range(0,2*len(tpos1)):
            res_tindices_list.remove(None) 
        
        res = _scalar(self.STbundle,res_tindices_list)
        return super()._trace_(res,tpos1,tpos2)

          
    def contract(self,*args): #contraction on both sindices or tindices
        
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg,tindices_object):
                other = arg
                it = i
                break
        else:
            raise TypeError("a tindices object must be provided in the argument list")
            
        if(it == 0):
            tpos1 = self._nid-1 #last as default
        elif(it == 1):   
            tpos1 = args[0]
        else:
            raise TypeError("Too many arguments")

        if(nargs-it-1 == 0):
            tpos2 = 0 #first as default
        elif(nargs-it-1 == 1):
            tpos2 = args[it+1]
        else:
            raise TypeError("Too many arguments")
        
 
        tpos1,tpos2 = self.check_tindices(other,tpos1,tpos2,"contract")
        res_tindices_list = list(self.tindices+other.tindices)
        
        for k in tpos1:
            res_tindices_list[k] = None
        for k in tpos2:
            res_tindices_list[k+self._nid] = None
        for k in range(0,2*len(tpos1)):
            res_tindices_list.remove(None)
         
        if(isinstance(other,_scalar)):
            res = _scalar(self.STbundle,res_tindices_list)
        elif(isinstance(other,_spin_tensor)):
            res = _spin_tensor(self.STbundle,res_tindices_list,other.sindices)
            
        return super()._contract_(other,res,tpos1,tpos2) # (*)

    
    def up_down(self,pos,u_d):
        
        if(self._nid == 0 ): raise TypeError("Object passed has no tindices")
        if(self.tindices[pos] == u_d): raise TypeError(f"Index in position {pos} is already {u_d}")

        eta = self.STbundle.eta(u_d)
        res = eta*self
        res = res.swap_tindices(0,pos+2,"move")
        res = res.trace(0,1)
        return res
    
    def up(self,pos=0):
        return self.up_down(pos,"up")
    
    def down(self,pos=0):
        return self.up_down(pos,"down")
    
    def par_der(self):
        res = self.STbundle.scalar(self.tindices+["down"])
        if(self._nid == 0):
            for k in range(0,4):
                res[k] = (self.comp.derivative())(self.STbundle.tframe[k])
            return res
        else:
            for ind, value in self._comp.items():
                for k in range(0,4):
                    res[ind + (k,)] = (value.derivative())(self.STbundle.tframe[k])
            return res

    def tcov_der(self):
        
        Ric_coef = self.STbundle.fb_coef
        
        res = self.par_der()

        for k in range(0,self._nid):
            if(self.tindices[k] == "up"):
                par = (self*Ric_coef).swap_tindices(k,self._nid,"move").trace(self._nid,self._nid+1)
                res += par
            elif(self.tindices[k] == "down"):
                par = (self*Ric_coef).swap_tindices(k,self._nid+1,"move").trace(self._nid,self._nid+1)
                res += par
        
        return res


from ST_bundle._spin_tensor import _spin_tensor

