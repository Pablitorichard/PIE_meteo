from .state import State

class History():

    def __init__(self, state_list):
        self.state_list = state_list
        self.size = len(state_list)

    @classmethod
    def fromCDF(cls, netCDF_file):
        size = netCDF_file.dimensions['Nt'].size 
        if (netCDF_file.Nt is not None):
            state_list = [State.fromCDF(netCDF_file, k) for k in range(size)]
            return cls(state_list)
        else:
            raise 'Empty CDF while initialising History'


    def save(self, netCDF_file, backup=False):
        for state in self.state_list:
            state.save(netCDF_file, backup)

    def append(self, state):
        self.state_list.append(state)
        self.size += 1

    def pop(self, k=None):
        try:
            self.state_list.pop(k)
            self.size -= 1
        except Exception as error:
            raise error

    #def forward(self, method):
    #    new_state = method(self)
    #    self.state_list.append(new_state)
    #    self.state_list.pop(0)