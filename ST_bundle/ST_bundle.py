from ST_bundle._scalar import _scalar
from ST_bundle._spin_tensor import _spin_tensor
from ST_bundle._sconnection import _sconnection
from ST_bundle._tconnection import _tconnection
from ST_bundle._totconnection import _totconnection
from sage.tensor.modules.comp import Components
from sage.rings.integer import Integer
from sage.functions.other import sqrt

class ST_bundle():
    
    def __init__(self,M,tframe,Orthonormal ="Yes"):

        self._M = M
        self._fbundle = M.vector_bundle(4, 'F', field ='complex')
        self._sbundle = M.vector_bundle(4, 'S', field='complex')
        self._tframe = [tframe,Orthonormal]
        self._tcoframe = tframe.coframe()
        self._fframe = self._fbundle.local_frame('e')
        self._sframe = self._sbundle.local_frame('e')
        self._metric = M.metric()
        self._invmetric = M.metric().inverse()
        self._connection = M.metric().connection()
        self._eta_up = self.eta_calculator("up")
        self._eta_down = self.eta_calculator("down")
        self._Gamma_up = self.Dirac_matrices_calculator()
        
        self._fb_coef = self.Ricci_rotation_coefficents(["up","down","down"])
        self._spinor_coef = self.spinor_coef_calculator()
    
    def __repr__(self):
        return f"ST bundle with respect to the Manifold {self.M} and to the tetrad frame {self.tframe}"

    @property
    def M(self):
        return self._M

    @property    
    def g(self):
        return self._metric
    
    @property
    def g_inv(self):
        return self._invmetric
    
    @property
    def connection(self):
        return self._connection
    
    @property
    def tframe(self):
        return self._tframe[0]
    
    @property
    def tcoframe(self):
        return self._tcoframe 
    
    @property
    def fframe(self):
        return self._fframe
    
    @property
    def sframe(self):
        return self._sframe
    
    @property
    def fbundle(self):
        return self._fbundle
    
    @property
    def sbundle(self):
        return self._sbundle

    @property   
    def fb_coef(self):
        return self._fb_coef
   
    @property
    def spinor_coef(self):
        return self._spinor_coef
   
    def scalar(self,tindices_list,scalar_field = None,state = "Mutable"):
        return _scalar(self,tindices_list,scalar_field,state)
    
    def spin_tensor(self,tindices_list,sindices_list,tensor = None,state = "Mutable"):
        return _spin_tensor(self,tindices_list,sindices_list,tensor,state)
    
    def tconnection(self):
        return _tconnection(self)
    
    def sconnection(self,A = None):
        return _sconnection(self,A)
    
    def totconnection(self,A = None):
        return _totconnection(self,A)

    
    def set_scalar_from_Components(self,tindices_list,Comp):
        
        if(Comp._frame != self.tframe and Comp._frame != self.fframe): raise TypeError("Components must be with respect to the tframe")
        
        if(not isinstance(Comp,Components)): raise TypeError("Required argument is not a Components object")
        index_list = list(Comp.index_generator())
        chart = Comp[index_list[0]].chart()
        if(len(index_list[0]) != len(tindices_list)): raise TypeError("Number of indices mismatch")
        
        res = self.scalar(tindices_list)
 
        for l in index_list:
            res[l] = self.M.scalar_field(Comp[l])
            
        res.set_immutable()
            
        return res
    
    
    def eta_calculator(self,index_type):
        
        Orthonormal = self._tframe[1]
        
        if(Orthonormal == "Yes"):
            eta = self.scalar([index_type,index_type])
            eta[0,0] = -self.M.scalar_field(1)
            eta[1,1] = self.M.scalar_field(1)
            eta[2,2] = self.M.scalar_field(1)
            eta[3,3] = self.M.scalar_field(1)
            
            eta.set_immutable()
            
            return eta
        
        elif(Orthonormal == "No"):
            eta = self.scalar([index_type,index_type])
            if(index_type == "down"):
                for i in range(0,4):
                    for j in range(0,4):
                        eta[i,j] = (self.g)(self.tframe[i],self.tframe[j])
            elif(index_type == "up"):
                for i in range(0,4):
                    for j in range(0,4):
                        eta[i,j] = (self.g_inv)(self.tcoframe[i],self.tcoframe[j])
            eta.set_immutable()           
                        
            return eta
    
    def eta(self,index_type):
        if(index_type == "up"): return self._eta_up
        elif(index_type == "down"): return self._eta_down
    
    def Ricci_rotation_coefficents(self,tindices_list):
        
        if(len(tindices_list) != 3): raise TypeError("Wrong number of indices")
        res = self.set_scalar_from_Components(["up","down","down"],self.connection.coef(self.tframe))
        
        if(tindices_list == ["down","down","up"]): res = res.down(0).up(2); res.set_immutable()
        if(tindices_list == ["down","up","down"]): res = res.down(0).up(1); res.set_immutable()
        
        if(tindices_list == ["up","up","down"]): res = res.up(1); res.set_immutable()
        if(tindices_list == ["down","up","up"]): res = res.down(0).up(1).up(2); res.set_immutable()
        if(tindices_list == ["up","down","up"]): res = res.up(2); res.set_immutable()
        
        if(tindices_list == ["up","up","up"]): res = res.up(1).up(2); res.set_immutable()
        if(tindices_list == ["down","down","down"]): res = res.down(0); res.set_immutable()
        
        return res
    
    
    def Zero(self,t_indices_list):
        res = self.spin_tensor(t_indices_list,(1,1),tensor = None,state = "Immutable")
        return res
    

    def Identity(self,t_indices_list):
        
        Id = self.sbundle.section_module().tensor((1,1))
        for k in range(0,4):
            Id[k,k] = 1
        res =  self.spin_tensor(t_indices_list,(1,1),Id)
        res.set_immutable()
        return res
    
  
    def Dirac_matrices_calculator(self):
        
        g0 = self.sbundle.section_module().tensor((1,1)) #The gamma matrices are defined as automorphism of the section module of the spinor bundle
        g1 = self.sbundle.section_module().tensor((1,1))
        g2 = self.sbundle.section_module().tensor((1,1)) 
        g3 = self.sbundle.section_module().tensor((1,1)) 

        g0[self.sframe,:] = ((1,0,0,0),(0,1,0,0),(0,0,-1,0),(0,0,0,-1)) #use I instead of i otherwise it is interpreted as an index
        g1[self.sframe,:] = ((0,0,0,1),(0,0,1,0),(0,-1,0,0),(-1,0,0,0)) # (+,-,-,-) signature
        g2[self.sframe,:] = ((0,0,0,-sqrt(-1)),(0,0,sqrt(-1),0),(0,sqrt(-1),0,0),(-sqrt(-1),0,0,0))
        g3[self.sframe,:] = ((0,0,1,0),(0,0,0,-1),(-1,0,0,0),(0,1,0,0))

        res = self.spin_tensor(["up"],(1,1))
        res[:] = [g0,g1,g2,g3]
        res.set_immutable()
        return res
        
    def Dirac_matrices(self,u_d):
        if(u_d == "up"):
            return self._Gamma_up
        elif(u_d == "down"):
            res = self._Gamma_up.down()
            res.set_immutable()
            return res

    def lgroup_spinor_generators(self,u_d):
        Gamma = self.Dirac_matrices(u_d)
        res = (Gamma.comm(Gamma))*(sqrt(-1)/2)
        if(u_d == "down"):
            res = res.down(0).down(1)
        res.set_immutable()
        return res
    
    def spinor_coef_calculator(self):
        
        Gamma = self.Dirac_matrices("up")
        Ric_coef_ddd = self.fb_coef.down(0)
        res = (Integer(1)/4)*Ric_coef_ddd.contract([0,1],Gamma@Gamma,[0,1]) #Gamma_{(a)}^{i}_{j} connection
        res.set_immutable()
        return res

