import abc
from sage.rings.integer import Integer

class tindices_object(abc.ABC):
    
    @staticmethod
    def recursive_generator(dim,tindices_list,gen):
        if len(tindices_list) == 0: return gen()
        return [tindices_object.recursive_generator(dim,tindices_list[1:],gen) for _ in range(dim)]

    @abc.abstractmethod
    def __init__(self,STb,tindices_list,object_type,state):

        self._STb = STb
        self._tindices = tindices_list
        self._comptype = type(object_type())
        self._comp = tindices_object.recursive_generator(4,tindices_list,object_type)
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
        if(len(self.tindices) != 0): 
            raise ValueError("Private attribute")
        return self._comp

    @comp.setter
    def comp(self,newcomp):
        if(len(self.tindices) != 0): 
            raise ValueError("Private attribute")	
        else:
            if(self.comp_type != type(newcomp)): raise TypeError("Wrong type assignment")
            if(self.state == "Immutable"): raise ValueError("the components of an immutable element cannot be changed")
            self._comp = newcomp
    
    @property
    def state(self):
        return self._state

    @state.setter
    def state(self,new_state):
        if(new_state != "Mutable" and new_state != "Immutable"): raise TypeError("State must be 'Mutable or 'Immutable'")
        if(self.state == "Immutable"): raise TypeError("Object is already immutable")
        self._state = new_state

    
    def set_immutable(self):
                
        self.state = "Immutable"
        lst,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        if(len(self.tindices) == 0):
            loop_command += f"self.comp.set_immutable()"
        else:
            loop_command += f"self[{lst[:-1]}].set_immutable()"            
        exec(loop_command)
        
            
    def __setitem__(self,pos_index,elem):
        
        if(self.state == "Immutable"): raise ValueError("the components of an immutable element cannot be changed")
        
        text = "self._comp"
        if(isinstance(pos_index,Integer) or isinstance(pos_index,int)):
            text += f"{[pos_index]}"
        else:
            for i in pos_index:
                text += f"[{i}]"
        text += " = elem"
        exec(text)
             
    def __getitem__(self,pos_index):
        
        text = "self._comp"
        if(isinstance(pos_index,Integer) or isinstance(pos_index,int)):
            text += f"{[pos_index]}"
        else:
            for i in pos_index:
                text += f"[{i}]"
        return eval(text)
        
    @staticmethod
    def nested_loop(index_name,start,stop,lst,loop_command,indent):
        
        for k in range(start,stop):
            lst += f"{index_name}{k},"
            loop_command += f"for {index_name}{k} in range(0,4):\n"
            indent += "\t"
            loop_command += indent
            
        return lst,loop_command,indent

    
    #Algorithms
    
    def copy_alg(self,res):

        lst,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        if(len(self.tindices) == 0):
            loop_command += f"res.comp = self.comp.copy()"
        else:
            loop_command += f"res[{lst[:-1]}] = self[{lst[:-1]}].copy()"            
        exec(loop_command)
        return res
    
    def eq_alg(self,other):
        
        if(self.tindices != other.tindices): return False
        count = [0]
        lst,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        loop_command += f"if(self[{lst[:-1]}] == other[{lst[:-1]}]): count[0] += 0\n"
        loop_command += indent
        loop_command += f"else: count[0] += 1"
        exec(loop_command)
        if(count[0] == 0): return True
        else: return False
        
    def sum_alg(self,other,res,operator):
        
        if(self.comp_type != other.comp_type): raise TypeError("Can't sum components of different types")
        if(self.tindices != other.tindices): raise TypeError("Can't sum two objects with different sets of indices")
        
        lst,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        
        if(operator == "+"): loop_command += f"res[{lst[:-1]}] = self[{lst[:-1]}] + other[{lst[:-1]}]"
        elif(operator == "-"): loop_command += f"res[{lst[:-1]}] = self[{lst[:-1]}] - other[{lst[:-1]}]"
        else: raise TypeError("Wrong operator type")
        exec(loop_command)
        return res

    def swap_tindices_alg(self,res,pos1,pos2):
        
        lst,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
            
        lst_swap = f"{lst[:3*pos1]}i{pos2},{lst[3*(pos1+1):3*pos2]}i{pos1},{lst[3*(pos2+1):]}"
        loop_command += f"res[{lst[:-1]}] = self[{lst_swap[:-1]}].copy()"
        
        exec(loop_command)
        
        return res
        
    def tensor_product_alg(self,other,res,spos1 = 0,spos2 = 0,operator = "*"):
        
        if(operator == "*"):
            op = "*"
            end = ""
        elif(operator == "@"):
            op = f".contract({spos1},"
            end = f",{spos2})"
        else:
            return TypeError("Wrong operator type")
        
        if(not isinstance(other,tindices_object)):
            lst_self,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
            loop_command += f"res[{lst_self[:-1]}] = self[{lst_self[:-1]}]{op}other{end}"
            exec(loop_command)
            return res
        
        if(self.STbundle.tframe != other.STbundle.tframe): raise TypeError("Objects have components with respect to different frames")
    
        if(len(self.tindices) == 0 and len(other.tindices) == 0):
            if(operator == "*"):
                res.comp = self.comp*other.comp
            else:
                res.comp = (self.comp).contract(spos1,other.comp,spos2)
            return res
    
        lst_self,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"","","")
        lst_other,loop_command,indent = tindices_object.nested_loop("j",0,len(other.tindices),"",loop_command,indent)
       
        lst_res = lst_self + lst_other
        
        if(len(self.tindices) == 0):
            loop_command += f"res[{lst_res[:-1]}] = self.comp{op}other[{lst_other[:-1]}]{end}"
        elif(len(other.tindices) == 0):
            loop_command += f"res[{lst_res[:-1]}] = self[{lst_self[:-1]}]{op}other.comp{end}"
        else:
            loop_command += f"res[{lst_res[:-1]}] = self[{lst_self[:-1]}]{op}other[{lst_other[:-1]}]{end}"
         
        exec(loop_command)
        return res
    
    def check_tindices(self,other,pos1,pos2,cont_type):

        if(not isinstance(pos1,list)): pos1 = [pos1]
        if(not isinstance(pos2,list)): pos2 = [pos2]

        if(cont_type == "trace"):
            if(len(self.tindices) <= 1): raise TypeError("Can't trace object with zero or one tindices")
            if(len(pos1)+len(pos2) > len(self.tindices)): raise TypeError("Too many indices")
        elif(cont_type == "contract"):
            if(len(self.tindices) == 0 or len(other.tindices) == 0): raise TypeError("Can't contract object with no tindices")
            if(len(pos1)>len(self.tindices) or len(pos2)>len(other.tindices)): raise TypeError("Too many indices")
            
        if(len(pos1) != len(pos2)): raise TypeError("The set of tindices to contract must be of the same size ")
        count = 0

        for k in range(0,len(pos1)):
            if(self.tindices[pos1[k]] == other.tindices[pos2[k]]): count += 1

        if(count != 0): raise TypeError("Can't trace two  up or two down indices ")

        return pos1,pos2
    
    
    def trace_alg(self,res,tpos1,tpos2,spos1,spos2,operator = "ttrace"): #pos1 and pos2 must be lists
        

        if(operator == "ttrace"):
            op = ""
        elif(operator == "strace"):
            op = ".trace(spos1,spos2)"
        else:
            raise TypeError("Wrong operator type")

        if(len(self.tindices) == 0 or len(self.tindices) == 1): raise TypeError("Can't trace an object with zero or one tindices")
        
        lst_run,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices)-2*len(pos1),"","","")
        _,loop_command,indent = tindices_object.nested_loop("s",0,len(pos1),"",loop_command,indent)  
        lst_total = [k for k in range(0,len(self.tindices))] # [0,1,2,...,number of tindices-1]
        
        lst = ""
        
        for k in range(0,len(self.tindices)):
            lst += f"i{k},"
          
        for k in range(0,len(pos1)):
            lst = lst.replace(f"i{pos1[k]},",f"s{k},")
            lst_total.remove(pos1[k])
        for k in range(0,len(pos2)):
            lst = lst.replace(f"i{pos2[k]},",f"s{k},")
            lst_total.remove(pos2[k])
        for k in range(0,len(lst_total)):
            lst = lst.replace(f"i{lst_total[k]},",f"i{k},")

        if(len(res.tindices) == 0):
            loop_command += f"res.comp += self[{lst[:-1]}]{op}"
        else:
            loop_command += f"res[{lst_run[:-1]}] += self[{lst[:-1]}]{op}"
        
        exec(loop_command)
        
        return res
   
    def contract_alg(self,other,res,tpos1,tpos2,spos1 = 0,spos2 = 0,operator = "*"):
    
        if(operator == "*"):
            op = "*"
            end = ""
        elif(operator == "@"):
            op = f".contract({spos1},"
            end = f",{spos2})"
        else:
            return TypeError("Wrong operator type")
    
        if(len(self.tindices) == 0 or len(other.tindices) == 0): raise TypeError("Can't contract an object with zero tindices")
        
        lst_run,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices)-len(tpos1),"","","")
        lst_run,loop_command,indent = tindices_object.nested_loop("j",0,len(other.tindices)-len(tpos2),lst_run,loop_command,indent)
        _,loop_command,indent = tindices_object.nested_loop("s",0,len(tpos1),"",loop_command,indent)  
        lst_total_self = [k for k in range(0,len(self.tindices))] # [0,1,2,...,number of tindices-1]
        lst_total_other = [k for k in range(0,len(other.tindices))]
        
        lst_self = ""
        lst_other = ""
        
        for k in range(0,len(self.tindices)):
            lst_self += f"i{k},"
        for k in range(0,len(other.tindices)):
            lst_other += f"j{k},"
          
        for k in range(0,len(tpos1)):
            lst_self = lst_self.replace(f"i{tpos1[k]},",f"s{k},")
            lst_total_self.remove(tpos1[k])
        for k in range(0,len(tpos2)):
            lst_other = lst_other.replace(f"j{tpos2[k]},",f"s{k},")
            lst_total_other.remove(tpos2[k])
        for k in range(0,len(lst_total_self)):
            lst_self = lst_self.replace(f"i{lst_total_self[k]},",f"i{k},")
        for k in range(0,len(lst_total_other)):
            lst_other = lst_other.replace(f"j{lst_total_other[k]},",f"j{k},")

        if(len(res.tindices) == 0):
            loop_command += f"res.comp += self[{lst_self[:-1]}]{op}other[{lst_other[:-1]}]{end}"
        else:
            loop_command += f"res[{lst_run[:-1]}] += self[{lst_self[:-1]}]{op}other[{lst_other[:-1]}]{end}"

        exec(loop_command)
        
        return res
    
    
    def get_tangent_tensor_alg(self,res,s_ind_sum):
            
        lst_res,loop_command,indent = tindices_object.nested_loop("j",0,s_ind_sum,"","","")
        lst_self,loop_command,indent = tindices_object.nested_loop("i",0,len(self.tindices),"",loop_command,indent)
        
        if(s_ind_sum == 0):
            loop_command += f"res[0] += self[{lst_self[:-1]}]*("
        else:
            loop_command += f"res{lst_res.replace(',',']').replace('j','[j')} += self[{lst_self[:-1]}][{lst_res[:-1]}]*("
        for k in range(0,len(self.tindices)):
            if(self.tindices[k] == "up"):
                loop_command += f"self.STbundle.tframe[i{k}]*"
            elif(self.tindices[k] == "down"):
                loop_command += f"self.STbundle.tcoframe[i{k}]*"
        loop_command = loop_command[:-1]
        loop_command += ")"
 
        exec(loop_command)
        
        return res
   
