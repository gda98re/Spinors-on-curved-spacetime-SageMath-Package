from ST_bundle.tindices_object import tindices_object

from sage.tensor.modules.tensor_free_module import FreeModuleTensor

class _spin_tensor(tindices_object):
    
    @staticmethod
    def NULL_tensor(STb,sindices_list):
        res = STb.sbundle.section_module().tensor(sindices_list)
        loop_command = "res[STb.sframe"
        for i in range(0,sindices_list[0]+sindices_list[1]):
            loop_command += ",0"
        loop_command += "] = 0"
        exec(loop_command)
        return res
    
    def __init__(self,STb,tindices_list,sindices_list,tensor = None,state = "Mutable"):
        
        if(tensor == None):
            _tensor = lambda : _spin_tensor.NULL_tensor(STb,sindices_list)
        else:
            _tensor = lambda : tensor.copy()

        if(_tensor().tensor_type() != sindices_list): raise TypeError("Tensor indices mismatch")       

        super().__init__(STb,tindices_list,_tensor,state)
        self._sindices = sindices_list
    
    def __setitem__(self,pos_index,elem):
        if(not isinstance(elem,self.comp_type)): raise TypeError("Can't assign something that is not a tensor on the spin bundle section module")
        super().__setitem__(pos_index,elem)
        
    def __repr__(self):
        return f"{self.tindices} Tetrad index family of Tensors {self.sindices} on the spin bundle section module with respect to the tetrad frame {self.STbundle.tframe}"
        
    def parent(self):
        return tindices_object
        
    @property
    def sindices(self):
        return self._sindices
    
    def display(self):
        
        res = super().recursive_generator(4,super().tindices,lambda : None)
        
        lst,loop_command,_ = super().nested_loop("i",0,len(self.tindices),"","","")
        if(len(self.tindices) == 0):
            return self.comp[:]
        else:
            loop_command += f"res{lst.replace(',',']').replace('i','[i')} = self[{lst[:-1]}][:]"
            exec(loop_command)
            return res

        
    def get_tangent_tensor(self):
        
        u = self.tindices.count("up")
        d = self.tindices.count("down")
        
        tensor_field = self.STbundle.M.tensor_field(u,d)
        tensor_field[self.STbundle.M.default_chart().frame(),[0 for k in range(0,len(self.tindices))]] = 0
        
        res = _spin_tensor.recursive_generator(4,["" for k in range(0,sum(self.sindices))],lambda : tensor_field)
        super().get_tangent_tensor_alg(res,sum(self.sindices))
        
        return res
    
    def copy(self):
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        return super().copy_alg(res)
    
            
    def __eq__(self,other):
        if(isinstance(other,_spin_tensor) == False): return False
        return super().eq_alg(other)
    
    def __add__(self,other):
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        super().sum_alg(other,res,"+")
        return res
    
    def __sub__(self,other):
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        super().sum_alg(other,res,"-")
        return res
   
    def swap_tindices(self,pos1 = 0,pos2 = 1):
        
        res_tindices_list = self.tindices[:pos1] + [self.tindices[pos2]] + self.tindices[pos1+1:pos2] + [self.tindices[pos1]] + self.tindices[pos2+1:]
        res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
        super().swap_tindices_alg(res,pos1,pos2)
        return res
    
    def swap_sindices(self,pos1 = 0,pos2 = 1):
        
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        lst_s,loop_command,indent = tindices_object.nested_loop("j",0,sum(self.sindices),"","","")
        lst_t,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"",loop_command,indent)
        lst_s_swap = f"{lst_s[:3*pos1]}j{pos2},{lst_s[3*(pos1+1):3*pos2]}j{pos1},{lst_s[3*(pos2+1):]}"
        loop_command += f"res[{lst_t[:-1]}][{lst_s[:-1]}] = (self[{lst_t[:-1]}].copy())[{lst_s_swap[:-1]}]"

        exec(loop_command)
        
        return res
        
    
    def __mul__(self,other): #tensor product full (all type of indices, both tetrad and spin)
        
        if(not isinstance(other,tindices_object)): 
            res_tindices_list = self.tindices
            if(isinstance(other,FreeModuleTensor)):
                res_sindices_list = (self.sindices[0]+other.tensor_type()[0],self.sindices[1]+other.tensor_type()[1])
            else:
                res_sindices_list = self.sindices
        else:
            res_tindices_list = self.tindices + other.tindices
            if(isinstance(other,_scalar)):
                res_sindices_list = self.sindices
            else:
                res_sindices_list = (self.sindices[0]+other.sindices[0],self.sindices[1]+other.sindices[1])
      
        res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
        super().tensor_product_alg(other,res,"*")
        return res
        
    def __rmul__(self,other):
        return self.__mul__(other)
        
    def __matmul__(self,other): #matrix multiplication between spin indices (contraction) and tensor product between tindices
        
        if(not isinstance(other,tindices_object)): 
            res_tindices_list = self.tindices
        else:
            res_tindices_list = self.tindices + other.tindices
            
        if(isinstance(other,_scalar)):
            res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
            super().tensor_product_alg(other,res,"*")
            
        elif(isinstance(other,_spin_tensor)):
            res_sindices_list = (self.sindices[0]+other.sindices[0]-1,self.sindices[1]+other.sindices[1]-1)
            res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
            super().tensor_product_alg(other,res,"@")
        elif(isinstance(other,FreeModuleTensor)):
            res_sindices_list = (self.sindices[0]+other.tensor_type()[0]-1,self.sindices[1]+other.tensor_type()[1]-1)
            res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
            super().tensor_product_alg(other,res,"@")
 
        return res
        
    def __rmatmul__(self,other):
        return self.__matmul__(other) 
    
    def comm(self,other): #buggy for more that two indices objects
        res = self@other - (other@self).swap_tindices(0,1)
        return res
    
    def anticomm(self,other): #buggy for more that two indices objects
        res = self@other + (other@self).swap_tindices(0,1)
        return res
        
    def ttrace(self,pos1=0,pos2=1):
                
        pos1,pos2 = self.check_tindices(self,pos1,pos2,"trace")
        res_tindices_list = list(self.tindices)
        
        for k in pos1:
            res_tindices_list[k] = None
        for k in pos2:
            res_tindices_list[k] = None
        for k in range(0,2*len(pos1)):
            res_tindices_list.remove(None) 
        
        res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
        super().trace_alg(res,pos1,pos2)
        if(len(res.tindices) == 0): 
            return res.comp
        else:
            return res
        
    def strace(self,pos1,pos2):
        res = _spin_tensor(self.STbundle,self.tindices,(self.sindices[0]-1,self.sindices[1]-1))
        lst_t,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        if(len(self.tindices) == 0):
            loop_command += f"res.comp = self.comp.trace(pos1,pos2)"
        else:
            loop_command += f"res[{lst_t[:-1]}] = self[{lst_t[:-1]}].trace(pos1,pos2)"
        exec(loop_command)
        return res
            
    def tcontract(self,pos1,other,pos2):
                
        pos1,pos2 = self.check_tindices(other,pos1,pos2,"contract")
        res_tindices_list = self.tindices+other.tindices
        
        for k in pos1:
            res_tindices_list[k] = None
        for k in pos2:
            res_tindices_list[k+len(self.tindices)] = None
        for k in range(0,2*len(pos1)):
            res_tindices_list.remove(None)
        
        if(not isinstance(other,tindices_object)): raise TypeError("Can't contract with objects that are not tindices_object instances") 
        
        if(isinstance(other,_scalar)):
            res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
            super().contract_alg(other,res,pos1,pos2,"*")
        else:
            res_sindices_list = (self.sindices[0]+other.sindices[0]-1,self.sindices[1]+other.sindices[1]-1)
            res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
            super().contract_alg(other,res,pos1,pos2,"@")
            
        if(len(res.tindices) == 0): 
            return res.comp
        else:
            return res
        
    def Dirac_adj(self):
    
        if(self.sindices != (1,0)): raise TypeError("Dirac adjoint works only on sections of the spinor bundle")

        res = _spin_tensor(self.STbundle,self.tindices,(0,1))

        lst,loop_command,indent = _spin_tensor.nested_loop("i",0,len(self.tindices),"","","")
        loop_command += "psi_bar = self.STbundle.section_module().tensor((0,1))\n"
        loop_command += indent
        loop_command += "for k in range(0,4):"
        loop_command += indent
        if(len(self.tindices) == 0):
            loop_command += "\t" + f"psi_bar[self.STbundle.sframe,k] = self.comp[k].expr().conjugate()\n"
        else:
            loop_command += "\t" + f"psi_bar[self.STbundle.sframe,k] = self[{lst[:-1]}][k].expr().conjugate()\n"
        loop_command += indent
        loop_command += f"res[{lst[:-1]}] = (psi_bar*(self.STbundle.Dirac_matrices('up')[0])).trace()"
        print(loop_command)
        exec(loop_command)

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
        
        res = self.STbundle.spin_tensor(self.tindices+["down"],self.sindices)
        lst_t,loop_command,indent = _spin_tensor.nested_loop("i",0,len(res.tindices),"","","")
        lst_s,loop_command,indent = _spin_tensor.nested_loop("j",0,sum(res.sindices),"",loop_command,indent)
        if(len(self.tindices) == 0):
            loop_command += f"res[{lst_t[:-1]}][{lst_s[:-1]}] = (self.STbundle.M.scalar_field(self.comp[{lst_s[:-1]}]).derivative())(self.STbundle.tframe[{lst_t[-3:-1]}])"
        else:
            loop_command += f"res[{lst_t[:-1]}][{lst_s[:-1]}] = (self.STbundle.M.scalar_field(self[{lst_t[:-4]}][{lst_s[:-1]}]).derivative())(self.STbundle.tframe[{lst_t[-3:-1]}])"

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
                res -= par
        
        return res
    
    
    def scov_der(self,s_conn = None):
        
        if(s_conn == None):
            s_conn = self.STbundle.spinor_coef
        
        res = self.par_der()
        
        for k in range(0,sum(self.sindices)):
            if(k < self.sindices[0]): #upper indices (first k in indices of a tensor of type (k,l))
                par = ((self*s_conn).swap_sindices(k,self.sindices[0])).strace(self.sindices[0],sum(self.sindices)+1)
                res += par
            elif(k >= self.sindices[0]): #lower indices (last l in indices of a tensor of type (k,l))
                par = ((self*s_conn).swap_sindices(k+1,sum(self.sindices)+1)).strace(self.sindices[0],sum(self.sindices)+1)
                res -= par
        
        return res
    
    
    def totcov_der(self,s_conn = None):
        
        Ric_coef = self.STbundle.fb_coef

        if(s_conn == None):
            s_conn = self.STbundle.spinor_coef

        res = self.par_der()
                
        for k in range(0,len(self.tindices)):
            if(self.tindices[k] == "up"):
                par = ((self*Ric_coef).swap_tindices(k,len(self.tindices))).ttrace(len(self.tindices),len(self.tindices)+1)
                res += par
            elif(self.tindices[k] == "down"):
                par = ((self*Ric_coef).swap_tindices(k,len(self.tindices)+1)).ttrace(len(self.tindices),len(self.tindices)+1)
                res -= par

        for k in range(0,sum(self.sindices)):
            if(k < self.sindices[0]): #upper indices (first k in indices of a tensor of type (k,l))
                par = ((self*s_conn).swap_sindices(k,self.sindices[0])).strace(self.sindices[0],sum(self.sindices)+1)
                res += par
            elif(k >= self.sindices[0]): #lower indices (last l in indices of a tensor of type (k,l))
                par = ((self*s_conn).swap_sindices(k+1,sum(self.sindices)+1)).strace(self.sindices[0],sum(self.sindices)+1)
                res -= par
        return res


from ST_bundle._scalar import _scalar
