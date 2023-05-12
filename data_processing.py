import numpy as np

class data_processing :
    def removeNan(array):
        nanind = np.isnan(array)
        return array[~nanind]
    def removeNan2(input, output):
        nanind = np.isnan(input)
        return output[~nanind]

