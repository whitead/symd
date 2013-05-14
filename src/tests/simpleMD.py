class SimpleMD:
    """SimpleMD Engine set-up and output parser"""    

    def __init__(self, integrator="vverlet", thermostat="bussi", force="lj"):
        name = "%s_%s_%s" % (force, integrator, thermostat)
        runParams = { 'steps':100,
                      'time_step':0.005,
                      'velocity_seed': 435423,
                      'anderson_nu':100,
                      'bussi_taut':5,
                      'thermostat_seed':54344,
                      'rcut':3,
                      'print_period':1,
                      'mass':

        
