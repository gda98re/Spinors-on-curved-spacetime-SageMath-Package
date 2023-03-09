import abc
from sage.rings.integer import Integer
from sage.tensor.modules.tensor_free_module import FreeModuleTensor
from sage.manifolds.scalarfield import ScalarField

class tindices_object(abc.ABC):
    
    @staticmethod
    def recursive_generator(dim,tindices_list,gen):
        if len(tindices_list) == 0: return gen()
        return [tindices_object.recursive_generator(dim,tindices_list[1:],gen) for _ in range(dim)]

    #@abc.abstractmethod
    def __init__(self,STb,tindices_list,object_type,state = "Mutable",start_index = 0):

        self._STb = STb
        self._tindices = tindices_list
        self._nid = len(tindices_list)
        self._dim = 4
        self._sindex = start_index #start index
        self._comptype = object_type #fa le veci di _ring
        if(self._nid == 0):
            self._comp = None #case of zero components
        else:
            self._comp = {} # the dictionary of components, with the index tuples
                        # as keys
        self._state = "Mutable"
        if(state == "Immutable"): self.set_immutable()
    
    @property   
    def STbundle(self):
        return self._STb
    
    @property
    def tindices(self):
        return self._tindices
    
    @property
    def comp_type(self):
        return self._comptype

    @property
    def comp(self):
        if(self._nid != 0): 
            raise ValueError("Private attribute")
        if(self._comp == None):
            return self.zero()
        else:
            return self._comp

    @comp.setter
    def comp(self,newcomp):
        if(self._nid != 0): 
            raise ValueError("Private attribute")
        else:
            if(self.state == "Immutable"): raise ValueError("the components of an immutable element cannot be changed")
            if(not isinstance(newcomp,self.comp_type)): raise ValueError("The assigned object has the wrong type")
            if(isinstance(newcomp,FreeModuleTensor)):
                if(newcomp.tensor_type() != self.sindices): raise TypeError("Tensor indices mismatch") 
            self._comp = newcomp
    
    @property
    def state(self):
        return self._state

    @state.setter
    def state(self,new_state):
        if(new_state != "Mutable" and new_state != "Immutable"): raise TypeError("State must be 'Mutable or 'Immutable'")
        if(self.state == "Immutable"): raise TypeError("Object is already immutable")
        self._state = new_state


    def check_tindices(self,other,pos1,pos2,cont_type):

        if(not isinstance(pos1,list)): pos1 = [pos1]
        if(not isinstance(pos2,list)): pos2 = [pos2]

        if(cont_type == "trace"):
            if(self._nid <= 1): raise TypeError("Can't trace object with zero or one tindices")
            if(len(pos1)+len(pos2) > self._nid): raise TypeError("Too many indices")
        elif(cont_type == "contract"):
            if(self._nid == 0 or len(other.tindices) == 0): raise TypeError("Can't contract object with no tindices")
            if(len(pos1)>self._nid or len(pos2)>len(other.tindices)): raise TypeError("Too many indices")
            
        if(len(pos1) != len(pos2)): raise TypeError("The set of tindices to contract must be of the same size ")
        count = 0

        for k in range(0,len(pos1)):
            if(self.tindices[pos1[k]] == other.tindices[pos2[k]]): count += 1

        if(count != 0): raise TypeError("Can't trace two  up or two down indices ")

        return pos1,pos2
    
    
    def set_immutable(self): #ok
        self.state = "Immutable"
        if(self._nid == 0):
            self.comp.set_immutable()
        else:
            for ind in self._comp.keys():
                self._comp[ind].set_immutable()        

    def _new_instance(self):
        return tindices_object(self.STbundle,self.tindices,self.comp_type,"Mutable")

    def _copy_(self,res):
        if(self._nid == 0):
            res.comp = self.comp.copy()
        else:
            for ind, val in self._comp.items():
                res._comp[ind] = val.copy()
        return res
    
    @staticmethod
    def NULL_tensor(STb,sindices_list): #ok
        res = STb.sbundle.section_module().tensor(sindices_list)
        command = "res[STb.sframe"
        for i in range(0,sindices_list[0]+sindices_list[1]):
            command += ",0"
        command += "] = 0"
        exec(command)
        return res
    
    def zero(self): #ok
        if(self.comp_type == FreeModuleTensor):
            return tindices_object.NULL_tensor(self.STbundle,self.sindices)
        elif(self.comp_type == ScalarField):
            return self.STbundle.M.scalar_field(0)
    
    def _del_zeros(self):
        # The zeros are first searched; they are deleted in a second stage, to
        # avoid changing the dictionary while it is read
        zeros = []
        for ind, value in self._comp.items():
            if value == 0:
                zeros.append(ind)
        for ind in zeros:
            del self._comp[ind]
              
    def _check_indices(self, indices):

        if isinstance(indices, (int, Integer)):
            ind = (indices,)
        else:
            ind = tuple(indices)
        if len(ind) != self._nid:
            raise ValueError(f"wrong number of indices: {self._nid} expected, while {len(ind)} are provided")
            
        si = self._sindex
        imax = self._dim - 1 + si
        for k in range(self._nid):
            i = ind[k]
            if i < si or i > imax:
                raise IndexError(f"index out of range: {i} not in [{si}, {imax}]")
        return ind
    
    def __getitem__(self, args):
        
        # Determining from the input the list of indices and the format
        if isinstance(args, (int, Integer, slice)):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        elif len(args) == self._nid:
            indices = args
        else:
            raise ValueError("too many indices provided")
            
        if isinstance(indices, slice):
            return self._get_list(indices)
        else:
            ind = self._check_indices(indices)
            if ind in self._comp:
                return self._comp[ind]
            else:  # if the value is not stored in self._comp, it is zero:
                return self.zero() #DA IMPLEmentARE!!!!!!1!!
    
    def _get_list(self, ind_slice):

        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            return [self[i] for i in range(start, stop)]
        
        if ind_slice.start is not None or ind_slice.stop is not None:
            raise NotImplementedError(f"function [start:stop] not implemented for components with {self._nid} indices")
                                      
        res = [self._gen_list([i]) for i in range(si, nsi)]
                                      
        if self._nid == 2:
            # 2-dim case: convert to matrix for a nicer output
            from sage.matrix.constructor import matrix
            from sage.structure.element import parent
            from sage.categories.rings import Rings
            if parent(resu[0][0]) in Rings():
                return matrix(res)
        return res

    def _gen_list(self, ind):

        if len(ind) == self._nid:
            return self[ind]
        else:
            si = self._sindex
            nsi = si + self._dim
            return [self._gen_list(ind + [i]) for i in range(si, nsi)]

    def __setitem__(self, args, value):
 
        if(self.state == "Immutable"): raise ValueError("the components of an immutable element cannot be changed")

        # Determining from the input the list of indices and the format
        if isinstance(args, (int, Integer, slice)):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        elif len(args) == self._nid:
            indices = args
        else:
            raise ValueError("too many indices provided")
                                      
        if isinstance(indices, slice):
            self._set_list(indices, value)
        else:
            ind = self._check_indices(indices)
            if(not isinstance(value,self.comp_type)): raise ValueError("The assigned object has the wrong type")
            if(isinstance(value,FreeModuleTensor)):
                if(value.tensor_type() != self.sindices): raise TypeError("Tensor indices mismatch") 
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if zero_value:
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done (zero components are not stored):
                if ind in self._comp:
                    del self._comp[ind]
            else:
                self._comp[ind] = value


    def _set_list(self, ind_slice, values):

        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            for i in range(start, stop):
                self[i] = values[i-start]
        else:
            if ind_slice.start is not None or ind_slice.stop is not None:
                raise NotImplementedError(f"function [start:stop] not implemented for components with {self._nid} indices")
            for i in range(si, nsi):
                self._set_value_list([i], values[i-si])

                
    def _set_value_list(self, ind, val):

        if len(ind) == self._nid:
            self[ind] = val
        else:
            si = self._sindex
            nsi = si + self._dim
            for i in range(si, nsi):
                self._set_value_list(ind + [i], val[i-si])
                    
    def display(self, symbol, latex_symbol=None, index_positions=None,
                index_labels=None, index_latex_labels=None,
                only_nonzero=True, only_nonredundant=False):
   
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        si = self._sindex
        nsi = si + self._dim
        if latex_symbol is None:
            latex_symbol = symbol
        if index_positions is None: #converts tindices list in "uddduudu..."
            string = str(self.tindices).replace("[","").replace("]","").replace("'","").replace(", ","")
            string = string.replace("p","").replace("own","")
            index_positions = string
        elif len(index_positions) != self._nid:
            raise ValueError("the argument 'index_positions' must contain " +
                             "{} characters".format(self._nid))
        if index_labels is None:
            index_labels = [str(i) for i in range(si, nsi)]
        elif len(index_labels) != self._dim:
            raise ValueError("the argument 'index_labels' must contain " +
                             "{} items".format(self._dim))
        # Index separator:
        max_len_symbols = max(len(s) for s in index_labels)
        if max_len_symbols == 1:
            sep = ''
        else:
            sep = ','
        if index_latex_labels is None:
            index_latex_labels = index_labels
        elif len(index_latex_labels) != self._dim:
            raise ValueError("the argument 'index_latex_labels' must " +
                             "contain {} items".format(self._dim))
        if only_nonredundant:
            generator = self.non_redundant_index_generator()
        else:
            generator = self.index_generator()
        rtxt = ''
        rlatex = r'\begin{array}{lcl}'
        for ind in generator:
            ind_arg = ind
            val = self[ind_arg]
            # Check whether the value is zero, preferably via the
            # fast method is_trivial_zero():
            if hasattr(val, 'is_trivial_zero'):
                zero_value = val.is_trivial_zero()
            else:
                zero_value = val == 0
            if not zero_value or not only_nonzero:
                indices = ''  # text indices
                d_indices = '' # LaTeX down indices
                u_indices = '' # LaTeX up indices
                previous = None  # position of previous index
                for k in range(self._nid):
                    i = ind[k] - si
                    if index_positions[k] == 'd':
                        if previous == 'd':
                            indices += sep + index_labels[i]
                        else:
                            indices += '_' + index_labels[i]
                        d_indices += r'\,' + "(" + index_latex_labels[i] + ")"
                        u_indices += r'\phantom{{\, ({})}}'.format(index_latex_labels[i])
                        previous = 'd'
                    else:
                        if previous == 'u':
                            indices += sep + index_labels[i]
                        else:
                            indices += '^' + index_labels[i]
                        d_indices += r'\phantom{{\, ({})}}'.format(index_latex_labels[i])
                        u_indices += r'\,' + "(" + index_latex_labels[i] + ")"
                        previous = 'u'
                rtxt += symbol + indices + ' = {} \n'.format(val)
                if(self.comp_type == FreeModuleTensor):
                    rlatex += (latex_symbol + r'_{' + d_indices + r'}^{'
                               + u_indices + r'} & = & ' + latex(val[:]) + r'\\')
                elif(self.comp_type == ScalarField):
                    rlatex += (latex_symbol + r'_{' + d_indices + r'}^{'
                               + u_indices + r'} & = & ' + latex(val.expr()) + r'\\')
                else:
                    rlatex += (latex_symbol + r'_{' + d_indices + r'}^{'
                               + u_indices + r'} & = & ' + latex(val) + r'\\')
        if rtxt == '':
            # no component has been displayed
            rlatex = ''
        else:
            # closing the display
            rtxt = rtxt[:-1]  # remove the last new line
            rlatex = rlatex[:-2] + r'\end{array}'
        return FormattedExpansion(rtxt, rlatex)                
                
    def is_zero(self):

        if not self._comp:
            return True

        #!# What follows could be skipped since _comp should not contain
        # any zero value
        # In other words, the full method should be
        #   return self.comp == {}
        for val in self._comp.values():
            if not (val == 0):
                return False
        return True
  
    def _swap_adjacent_tindices_(self, res, pos1, pos2, pos3,ref):
        
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            if(ref == "move"):
                res._comp[new_ind] = val
            elif(ref == "copy"):
                res._comp[new_ind] = val.copy()
            else:
                raise ValueError("Wrong ref type")
            # the above writing is more efficient than res[new_ind] = val
            # it does not work for the derived class CompWithSym, but for the
            # latter, the function CompWithSym.swap_adjacent_indices will be
            # called and not the present function.
        return res
    
    def _swap_tindices_(self, res, pos1, pos2,ref):
        
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + (ind[pos2],) + ind[pos1+1:pos2] + (ind[pos1],) + ind[pos2+1:]
            if(ref == "move"):
                res._comp[new_ind] = val
            elif(ref == "copy"):
                res._comp[new_ind] = val.copy()
            else:
                raise ValueError("Wrong ref type")
            # the above writing is more efficient than res[new_ind] = val
            # it does not work for the derived class CompWithSym, but for the
            # latter, the function CompWithSym.swap_adjacent_indices will be
            # called and not the present function.
        return res

    def __eq__(self, other):

        if isinstance(other, (int, Integer)): # other is 0
            if other == 0:
                return self.is_zero()
            else:
                raise TypeError("cannot compare a set of components to a number")
        else: # other is another Components
            if not isinstance(other, tindices_object):
                raise TypeError("an instance of tindices_object is expected")
            if other.comp_type != self.comp_type:
                return False
            if other.STbundle._tframe != self.STbundle._tframe:
                return False
            if other._nid != self._nid:
                return False
            if other._sindex != self._sindex:
                return False
            return (self - other).is_zero()

    def __ne__(self, other): #ok
        return not self == other

    def _neg_(self,res): #- ok
        for ind, val in self._comp.items():
             res._comp[ind] = - val
        return res
    
    def _add_(self, other, res): #ok

        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, tindices_object):
            raise TypeError("the second argument for the addition must be " +
                            "an instance of tindices_object")
        #if isinstance(other, tindices_objectWithSym):
        #    return other + self     # to deal properly with symmetries
        if other.STbundle.tframe != self.STbundle.tframe:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._nid != self._nid:
            raise ValueError("the two sets of components do not have the " +
                             "same number of indices")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        # Initialization of the res to self.copy(), so that there remains
        # only to add other:
       
        # Sequential computation
        for ind, val in other._comp.items():
            res[ind] += val

        return res

    def __radd__(self, other): #ok
        return self + other

    def __sub__(self, other): #ok
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        return self + (-other)  #!# correct, deals properly with
                                # symmetries, but is probably not optimal
    def __rsub__(self, other): #ok
        return (-self) + other
    
    def _truediv_(self, other,res): #ok
        if isinstance(other, tindices_object):
            raise NotImplementedError("division by an object of type " +
                                      "tindices_object not implemented")
        for ind, val in self._comp.items():
            res._comp[ind] = val / other
        return res
    
    def _tensor_product_(self,other,res,spos1 = 0,spos2 = 0,operator="*"):
        if(not isinstance(other,tindices_object)):
            # Sequential computation
            if(self._nid == 0):
                if(operator == "*"):
                    res.comp = self.comp*other
                elif(operator == "@"):
                    res.comp = (self.comp).contract(spos1,other,spos2)
                else:
                    return TypeError("Wrong operator type")
                return res
            else:
                if(operator == "*"):
                    for ind_s, val_s in self._comp.items():
                        res._comp[ind_s] = val_s * other
                elif(operator == "@"):
                    for ind_s, val_s in self._comp.items():
                        res._comp[ind_s] = val_s.contract(spos1,other,spos2)
                else:
                    return TypeError("Wrong operator type")
                return res
        
        if(self.STbundle.tframe != other.STbundle.tframe): raise TypeError("Objects have components with respect to different frames")
    
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        
        # Sequential computation
        if(self._nid == 0 and other._nid == 0):
            if(operator == "*"):
                res.comp = self.comp*other.comp
            elif(operator == "@"):
                res.comp = (self.comp).contract(spos1,other.comp,spos2)
            else:
                return TypeError("Wrong operator type")
            return res
        elif(self._nid == 0):
            for ind_o, val_o in other._comp.items():
                if(operator == "*"):
                    res._comp[ind_o] = self.comp * val_o
                elif(operator == "@"):
                    res._comp[ind_o] = (self.comp).contract(spos1,val_o,spos2)
                else:
                    return TypeError("Wrong operator type")
            return res
        elif(other._nid == 0):
            for ind_s, val_s in self._comp.items():
                if(operator == "*"):
                    res._comp[ind_s] = val_s * other.comp
                elif(operator == "@"):
                    res._comp[ind_s] = val_s.contract(spos1,other.comp,spos2)
                else:
                        return TypeError("Wrong operator type")
            return res
        else:
            for ind_s, val_s in self._comp.items():
                for ind_o, val_o in other._comp.items():
                    if(operator == "*"):
                        res._comp[ind_s + ind_o] = val_s * val_o
                    elif(operator == "@"):
                        res._comp[ind_s + ind_o] = val_s.contract(spos1,val_o,spos2)
                    else:
                        return TypeError("Wrong operator type")
            return res
    
    def _trace_(self,res,tpos1,tpos2,spos1=0,spos2=0,operator = "ttrace"):

        si = self._sindex
        nsi = si + self._dim
        
        ncontr = len(tpos1) # number of contractions
        res_nid = self._nid - 2*ncontr
                                      
        if self._nid == 2:
            res = 0
            for i in range(si, nsi):
                if(operator == "ttrace"):
                    res += self[i,i]
                elif(operator == "strace"):
                    res += self[i,i].trace(spos1,spos2)
                else:
                    raise TypeError("Wrong operator type")
            return res
        
        else:
            pos_s = [None for i in range(self._nid)]  # initialization
            contractions = [(tpos1[i], tpos2[i]) for i in range(ncontr)]
            
            comp_for_contr = tindices_object(self.STbundle,["" for i in range(ncontr)],
                                 self.comp_type,state="Mutable",start_index=self._sindex)
            shift_o = self._nid - ncontr
            shift = 0
            for pos in range(self._nid):
                for contract_index in (tpos1+tpos2):
                    if pos == contract_index:
                        shift += 1
                        break
                else:
                        pos_s[pos] = pos - shift
            rev_s = [pos_s.index(i) for i in range(res._nid)]
            # sequential
            for ind in res.index_generator():
                ind_s = [None for i in range(self._nid)]  # initialization
                for i, pos in enumerate(rev_s):
                    ind_s[pos] = ind[i]
                sm = 0
                for ind_c in comp_for_contr.index_generator():
                    ic = 0
                    for pos1, pos2 in contractions:
                        k = ind_c[ic]
                        ind_s[pos1] = k
                        ind_s[pos2] = k
                        ic += 1
                    if(operator == "ttrace"):
                        sm += self[ind_s]
                    elif(operator == "strace"):
                        sm += self[ind_s].trace(spos1, spos2)
                    else:
                        return TypeError("Wrong operator type")
                res[ind] = sm

            return res

    def _contract_(self,other,res, tpos1, tpos2, spos1 = 0, spos2 = 0, operator = "*"):

        ncontr = len(tpos1) # number of contractions
        contractions = [(tpos1[i], tpos2[i]) for i in range(ncontr)]
        res_nid = self._nid + other._nid - 2*ncontr
        #
        # Special case of a scalar res
        #
        if res_nid == 0:
            # To generate the indices tuples (of size ncontr) involved in the
            # the contraction, we create an empty instance of Components with
            # ncontr indices and call the method index_generator() on it:
            comp_for_contr = tindices_object(self.STbundle,["" for i in range(ncontr)],
                                 self.comp_type,state="Mutable",start_index=self._sindex)
            # sequential
            res = 0
            for ind in comp_for_contr.index_generator():
                if(operator == "*"):
                    res += self[ind] * other[ind]
                elif(operator == "@"):
                    res += self[ind].contract(spos1, other[ind], spos2)
                else:
                    raise TypeError("Wrong operator type")
            return res
                                      
        # Positions of self and other indices in the res
        #  (None = the position is involved in a contraction and therefore
        #   does not appear in the final res)
        #
        comp_for_contr = tindices_object(self.STbundle,["" for i in range(ncontr)],
                                 self.comp_type,state="Mutable",start_index=self._sindex)
        shift_o = self._nid - ncontr
        pos_s = [None for i in range(self._nid)]  # initialization
        pos_o = [None for i in range(other._nid)] # initialization
        shift = 0
        for pos in range(self._nid):
            for contract_pair in contractions:
                if pos == contract_pair[0]:
                    shift += 1
                    break
            else:
                pos_s[pos] = pos - shift
        for pos in range(other._nid):
            for contract_pair in contractions:
                if pos == contract_pair[1]:
                    shift += 1
                    break
            else:
                    pos_o[pos] = self._nid + pos - shift
        rev_s = [pos_s.index(i) for i in range(self._nid-ncontr)] #rev_s is a list of the positions of the non contracted indices of self
        rev_o = [pos_o.index(i) for i in range(self._nid-ncontr, res_nid)] #rev_o is a list of the positions of the non contracted indices of other

        # sequential
        for ind in res.index_generator():
            ind_s = [None for i in range(self._nid)]  # initialization
            ind_o = [None for i in range(other._nid)] # initialization
            for i, pos in enumerate(rev_s):
                ind_s[pos] = ind[i]
            for i, pos in enumerate(rev_o):
                ind_o[pos] = ind[shift_o+i]
            sm = 0
            for ind_c in comp_for_contr.index_generator():
                ic = 0
                for pos_s, pos_o in contractions:
                    k = ind_c[ic]
                    ind_s[pos_s] = k
                    ind_o[pos_o] = k
                    ic += 1
                if(operator == "*"):
                    sm += self[ind_s] * other[ind_o]
                elif(operator == "@"):
                    sm += self[ind_s].contract(spos1,other[ind_o],spos2)
                else:
                    return TypeError("Wrong operator type")
            res[ind] = sm

        return res

    
    def get_tangent_tensor(self):
 
        u = self.tindices.count("up")
        d = self.tindices.count("down")
        NULL_tg_tensor = self.STbundle.M.tensor_field(u,d)
        NULL_tg_tensor[self.STbundle.M.default_chart().frame(),[0 for k in range(0,self._nid)]] = 0
        res = {}

        if(isinstance(self[next(self.index_generator())],ScalarField)):
            res[0] = NULL_tg_tensor
            command = "for ind,value in self._comp.items():\n"
            command += "\t"
            command += f"res[0] += value*("
            for k in range(0,self._nid):
                if(self.tindices[k] == "up"):
                    command += f"self.STbundle.tframe[ind[{k}]]*"
                elif(self.tindices[k] == "down"):
                    command += f"self.STbundle.tcoframe[ind[{k}]]*"
            command = command[:-1]
            command += ")"
            exec(command)
            return res[0]

        else:    
            for inds in self[next(self.index_generator())].comp().index_generator(): #init of the res
                res[inds] = NULL_tg_tensor.copy() 

            command = "for inds in self[next(self.index_generator())].comp().index_generator():\n"
            command += "\t"
            command += "for ind,value in self._comp.items():\n"
            command += "\t\t"
            command += f"res[inds] += value[inds]*("
            for k in range(0,self._nid):
                if(self.tindices[k] == "up"):
                    command += f"self.STbundle.tframe[ind[{k}]]*"
                elif(self.tindices[k] == "down"):
                    command += f"self.STbundle.tcoframe[ind[{k}]]*"
            command = command[:-1]
            command += ")"
            exec(command)
            return res

        
    def index_generator(self):
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(self._nid)]
        ind_end = [si for k in range(self._nid)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(self._nid-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def non_redundant_index_generator(self):
        for ind in self.index_generator():
            yield ind
   
