import src.results.util as ut
from joblib import load 


def get_target_crm(input,env,options,args):
    """
    Returns the set of predicted co-regulatory modules for the target in an experiment as a dataframe 
    """
    result_type = options.get("result_name","default")
    out_path, temp_path = ut.get_exp_path(input,env)




"""
To interactively explore an experiment
"""
class Explore:
    def __init__(self,env_path,exp_file_path):
        self.env = ut.get_env(env_path)
        self.input = ut.get_input(exp_file_path)
        output_path, temp_path  = ut.get_exp_path(self.input,self.env)
        self.output_path = output_path
        self.temp_path = temp_path
        self.default_args = {"interactive":True}
        self.target_cache = {}
        print("done")
    def _target_cache(self,target):
        if target in self.target_cache:
            return self.target_cache[target]
        else:
            tfile_path = self.temp_path / "results"/ f"{target}.pkl"
            tfile = load(tfile_path)
            self.target_cache[target] = tfile
    def target_gf(self,target):
        """"""
        data = self._target_cache(target)
        return data["result"]["gene_frequency"]

    def target_crm(self,target,result_name="default"):
        """"""
        options = { "result_type":result_name, "target":target }
        args = {"interactive":True}
        return get_target_crm(self.input,self.env,options,self.default_args)
    
    

