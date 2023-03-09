from ST_bundle.tindices_object import tindices_object

from sage.tensor.modules.tensor_free_module import FreeModuleTensor

class _spin_tensor(tindices_object):
    
    def __init__(self,STb,tindices_list,sindices_list,tensor = None, state = "Mutable"):
        self._sindices = sindices_list
        super().__init__(STb,tindices_list,FreeModuleTensor,state)
        if(tensor != None):
            _tensor = lambda : tensor.copy()
            lis = _spin_tensor.recursive_generator(4,tindices_list,_tensor)
            if(self._nid == 0):
                self.comp = lis
            else:
                self[:] = lis     

        
    def __repr__(self):
        return f"{self.tindices} Tetrad index family of Tensors {self.sindices} on the spin bundle section module with respect to the tetrad frame {self.STbundle.tframe}"
        
    def parent(self):
        return tindices_object
        
    @property
    def sindices(self):
        return self._sindices
    
    def copy(self):
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        return super()._copy_(res)
    
    def swap_adjacent_tindices(self, pos1, pos2, pos3, ref = "copy"):
        res_tindices_list = self.tindices[:pos1] + self.tindices[pos2:pos3] + self.tindices[pos1:pos2] + self.tindices[pos3:]
        res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
        return super()._swap_adjacent_tindices_(res,pos1, pos2, pos3,ref)
    
    def swap_tindices(self, pos1 = 0, pos2 = 1, ref = "copy"):
        res_tindices_list = self.tindices[:pos1] + [self.tindices[pos2]] + self.tindices[pos1+1:pos2] + [self.tindices[pos1]] + self.tindices[pos2+1:]
        res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
        return super()._swap_tindices_(res,pos1,pos2,ref)
    
    def swap_sindices(self,pos1 = 0,pos2 = 1, ref = "copy"):
        
        res = _spin_tensor(self.STbundle,self.tindices,self.sindices)
        for ind, value in self._comp.items():
            par = self.zero()
            for inds in self[next(self.index_generator())].comp().index_generator():
                new_inds = inds[:pos1] + (inds[pos2],) + inds[pos1+1:pos2] + (inds[pos1],) + inds[pos2+1:]
                if(ref == "copy"):
                    par[new_inds] = value[inds].copy()
                elif(ref == "move"):
                    par[new_inds] = value[inds]
                else:
                    raise ValueError("Wrong ref type")
            res[ind] = par
        return res
    
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
        return super()._tensor_product_(other,res) # (*)

    def __rmul__(self,other):
        return self.__mul__(other)

    def __matmul__(self,other): #matrix multiplication between spin indices (contraction) and tensor product between tindices
        
        if(not isinstance(other,tindices_object)): 
            res_tindices_list = self.tindices
        else:
            res_tindices_list = self.tindices + other.tindices
            
        if(isinstance(other,_scalar)):
            res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
            return super()._tensor_product_(other,res) # (*)
            
        elif(isinstance(other,_spin_tensor)):
            res_sindices_list = (self.sindices[0]+other.sindices[0]-1,self.sindices[1]+other.sindices[1]-1)
            res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
            return super()._tensor_product_(other,res,sum(self.sindices)-1,0,"@")
        elif(isinstance(other,FreeModuleTensor)):
            res_sindices_list = (self.sindices[0]+other.tensor_type()[0]-1,self.sindices[1]+other.tensor_type()[1]-1)
            res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
            return super()._tensor_product_(other,res,sum(self.sindices)-1,0,"@")
        else:
            raise TypeError("Wrong type")

    def __rmatmul__(self,other):
        return self.__mul__(other)

    def __truediv__(self,other):
        res = self.copy()
        return super()._truediv_(other,res)
    
    
    def comm(self,other):
        res = self@other - (other@self).swap_adjacent_tindices(0,self._nid,self._nid+other._nid,"move")
        return res
    
    def anticomm(self,other):
        res = self@other + (other@self).swap_adjacent_tindices(0,self._nid,self._nid+other._nid,"move")
        return res
    
    def trace(self,*args): #trace on both sindices or tindices
        
        args = list(args) #I want to use remove
        nargs = len(args)
        if(nargs == 0): 
            if(self.tindices == ["up","down"] or self.tindices == ["down","up"]):
                tpos1 = 0
                tpos2 = 1
                spos1 = None
                spos2 = None
            else:
                raise TypeError("At least two indices must be provided")
        elif(nargs == 1):
            if(args[0] == None):
                tpos1 = None
                tpos2 = None
                spos1 = self.sindices[0]-1 #last of the upper sindices as default
                spos2 = self.sindices[0] #first of the lower sindices as default
            else:
                raise TypeError("At least two indices must be provided")
        elif(nargs == 2):
            tpos1 = args[0]
            tpos2 = args[1]
            spos1 = None
            spos2 = None
        elif(nargs == 3):
            if (None in args):
                args.remove(None)
                spos1 = args[0]
                spos2 = args[1]
                tpos1 = None
                tpos2 = None
            else: raise TypeError("Can't pass three indices as arguments")
        elif(nargs == 4):
            tpos1 = args[0]
            tpos2 = args[1]
            spos1 = args[2]
            spos2 = args[3]
        else:
            raise TypeError("Too many arguments")    
            
        res_tindices_list = self.tindices
        
        if(tpos1 != None and tpos2 != None):   #####tindices trace
            
            tpos1,tpos2 = self.check_tindices(self,tpos1,tpos2,"trace")
            #checks on spos1 and spos2 are done by trace() method in FreeTensorModule class
            res_tindices_list = list(self.tindices)

            for k in tpos1:
                res_tindices_list[k] = None
            for k in tpos2:
                res_tindices_list[k] = None
            for k in range(0,2*len(tpos1)):
                res_tindices_list.remove(None) 

            if((spos1 == None or spos2 == None)): #####tindices trace only
                res = _spin_tensor(self.STbundle,res_tindices_list,self.sindices)
                return super()._trace_(res,tpos1,tpos2) # (*)
            else:                                 #####trace on both type of indices
                res_sindices_list = (self.sindices[0]-1,self.sindices[1]-1)

                if(res_sindices_list == (0,0)): #the contraction of sindices produces a scalar field
                        res = _scalar(self.STbundle,res_tindices_list)
                else:
                        res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
                return super()._trace_(res,tpos1,tpos2,spos1,spos2,"strace")
        
        
        else: #if tpos1 or tpos2 are None -> no trace on tindices
            
            if(spos1 != None and spos2 != None): #sindices trace only
                
                if((self.sindices[0]-1,self.sindices[1]-1) == (0,0)): #the contraction of sindices produces a scalar field
                        res = _scalar(self.STbundle,self.tindices)
                else:
                        res = _spin_tensor(self.STbundle,self.tindices,(self.sindices[0]-1,self.sindices[1]-1))
                if(self._nid == 0):
                    res.comp = self.comp.trace(spos1,spos2)
                    return res.comp
                else:
                    for ind, value in self._comp.items():
                        res[ind] = value.trace(spos1,spos2)
                    return res
            else:
                raise TypeError("Must provide tindices to trace")
            
        
    def contract(self,*args): #contraction on both sindices or tindices
        
        nargs = len(args)
        for i, arg in enumerate(args):
            if(isinstance(arg,tindices_object)):
                other = arg
                it = i
                break
            elif(isinstance(arg,FreeModuleTensor)):
                other = _spin_tensor(self.STbundle,[],arg.tensor_type(),arg)
                it = i
                break
        else:
            raise TypeError("a tindices object must be provided in the argument list")
        if(it == 0):
            tpos1 = self._nid-1 #last as default
            spos1 = sum(self.sindices)-1 #last as default
        elif(it == 1):   
            tpos1 = args[0]
            spos1 = sum(self.sindices)-1
        elif(it == 2):
            if(isinstance(other,_scalar)): raise TypeError("sindices provided, but other is a scalar object")
            tpos1 = args[0]
            spos1 = args[1]
        else:
            raise TypeError("Too many arguments")

        if(nargs-it-1 == 0):
            tpos2 = 0 #first as default
            spos2 = 0 #first as default
        elif(nargs-it-1 == 1):
            tpos2 = args[it+1]
            spos2 = 0
        elif(nargs-it-1 == 2):
            if(isinstance(other,_scalar)): raise TypeError("sindices provided, but other is a scalar object")
            tpos2 = args[it+1]
            spos2 = args[it+2]
        else:
            raise TypeError("Too many arguments")
            
        if(isinstance(other,_scalar)): spos1,spos2 = None, None    
            
        res_tindices_list = self.tindices+other.tindices
        
        
        if(tpos1 != None and tpos2 != None):   #####tindices contraction
            
            tpos1,tpos2 = self.check_tindices(other,tpos1,tpos2,"contract")
            #checks on spos1 and spos2 are done by contract() method in FreeTensorModule class
            for k in tpos1:
                res_tindices_list[k] = None
            for k in tpos2:
                res_tindices_list[k+self._nid] = None
            for k in range(0,2*len(tpos1)):
                res_tindices_list.remove(None)
        
            if((spos1 == None or spos2 == None)): #####tindices contraction only
                
                if(isinstance(other,_scalar)):
                    res_sindices_list = self.sindices
                else:
                    res_sindices_list = (self.sindices[0]+other.sindices[0],self.sindices[1]+other.sindices[1])
                    
                res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)
                return super()._contract_(other,res,tpos1,tpos2) # (*)
                
            else:                #####contraction on both type of indices
                
                res_sindices_list = (self.sindices[0]+other.sindices[0]-1,self.sindices[1]+other.sindices[1]-1)

                if(res_sindices_list == (0,0)):
                    res = _scalar(self.STbundle,res_tindices_list)
                else:
                    res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)

                return super()._contract_(other,res,tpos1,tpos2,spos1,spos2,"@")
        
        
        else: #if tpos1 or tpos2 are None -> no contraction on tindices
            
            if(spos1 != None and spos2 != None): #sindices contraction only
            
                res_sindices_list = (self.sindices[0]+other.sindices[0]-1,self.sindices[1]+other.sindices[1]-1)

                if(res_sindices_list == (0,0)):
                    res = _scalar(self.STbundle,res_tindices_list)
                else:
                    res = _spin_tensor(self.STbundle,res_tindices_list,res_sindices_list)

                return super()._tensor_product_(other,res,spos1,spos2,"@")
            
            else:
                raise TypeError("Must provide tindices to contract")
                  

        
    def Dirac_adj(self):
    
        if(self.sindices != (1,0)): raise TypeError("Dirac adjoint works only on sections of the spinor bundle")
        res = _spin_tensor(self.STbundle,self.tindices,(0,1))
        psi_bar = self.STbundle.sbundle.section_module().tensor((0,1))
        if(self._nid == 0): 
            for k in range(0,4):
                psi_bar[self.STbundle.sframe,k] = self.comp[k].expr().conjugate()
                res.comp = (psi_bar*(self.STbundle.Dirac_matrices('up')[0])).trace()
            return res

        for ind, value in self._comp.items():
            for k in range(0,4):
                psi_bar[self.STbundle.sframe,k] = value[k].expr().conjugate()
                res[ind] = (psi_bar*(self.STbundle.Dirac_matrices('up')[0])).trace()
        return res
       

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
        
        res = self.STbundle.spin_tensor(self.tindices+["down"],self.sindices)
        if(self._nid == 0):
            for k in range(0,4):
                par = self.zero()
                for inds in self.comp.comp().index_generator():
                    par[inds] = (self.STbundle.M.scalar_field(self.comp[inds]).derivative())(self.STbundle.tframe[k])
                res[k] = par
            return res
        else:
            for ind, value in self._comp.items():
                for k in range(0,4):
                    par = self.zero()
                    for inds in self[next(self.index_generator())].comp().index_generator():
                        par[inds] = (self.STbundle.M.scalar_field(value[inds]).derivative())(self.STbundle.tframe[k])
                    res[ind+(k,)] = par
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
                res -= par
        
        return res
    
    
    def scov_der(self,s_conn = 0):
        
        if(s_conn == 0):
            s_conn = self.STbundle.spinor_coef
        
        res = self.par_der()

        #if(self.sindices == (1,1)):
        #    res += -self.comm(s_conn)
        #    return res
        
        for k in range(0,sum(self.sindices)):
            if(k < self.sindices[0]): #upper indices (first k in indices of a tensor of type (k,l))
                par = (self*s_conn).swap_sindices(k,self.sindices[0],"move").trace(None,self.sindices[0],sum(self.sindices)+1)
                res += par
            elif(k >= self.sindices[0]): #lower indices (last l in indices of a tensor of type (k,l))
                par = (self*s_conn).swap_sindices(k+1,sum(self.sindices+1),"move").trace(None,self.sindices[0],sum(self.sindices)+1)
                res -= par
        
        return res
    
    
    def totcov_der(self,s_conn = 0):
        
        Ric_coef = self.STbundle.fb_coef

        if(s_conn == 0):
            s_conn = self.STbundle.spinor_coef

        res = self.par_der()
                
        for k in range(0,self._nid):
            if(self.tindices[k] == "up"):
                par = (self*Ric_coef).swap_tindices(k,self._nid,"move").trace(self._nid,self._nid+1)
                res += par
            elif(self.tindices[k] == "down"):
                par = (self*Ric_coef).swap_tindices(k,self._nid+1,"move").trace(self._nid,self._nid+1)
                res -= par

        if(self.sindices == (1,1)):
            res += -self.comm(s_conn)
            return res

        for k in range(0,sum(self.sindices)):
            if(k < self.sindices[0]): #upper indices (first k in indices of a tensor of type (k,l))
                par = (self*s_conn).swap_sindices(k,self.sindices[0],"move").trace(None,self.sindices[0],sum(self.sindices)+1)
                res += par
            elif(k >= self.sindices[0]): #lower indices (last l in indices of a tensor of type (k,l))
                par = (self*s_conn).swap_sindices(k+1,sum(self.sindices)+1,"move").trace(None,self.sindices[0],sum(self.sindices)+1)
                res -= par
        return res


from ST_bundle._scalar import _scalar
