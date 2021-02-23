from .state import State

class History():

    def __init__(self, state_list):
        self.state_list = state_list
        self.size = len(state_list)

    @classmethod
    def fromCDF(cls, netCDF_file):
        size = netCDF_file.dimensions['Nt'].size 
        if size > 0:
            state_list = [State.fromCDF(netCDF_file, k) for k in range(size)]
            return cls(state_list)
        else:
            raise 'Empty CDF while initialising History'


    def save(self, netCDF_file, backup=False, saved_variables=None):
        if not backup:
        # the oldest state is supposed to be completely known
        # thus it is this one we choose to export
            self.state_list[0].save(netCDF_file, backup=False, saved_vrs=saved_variables)
        else:
            for ind, state in enumerate(self.state_list):
                state.save(netCDF_file, backup=True, k=ind)

    def append(self, state):
        self.state_list.append(state)
        self.size += 1

    def pop(self, k=None):
        try:
            self.state_list.pop(k)
            self.size -= 1
        except Exception as error:
            raise error
