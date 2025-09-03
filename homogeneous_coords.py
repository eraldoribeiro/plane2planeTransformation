# homogeneous_coords.py
import numpy as np

class HomogeneousCoords:

    @staticmethod
    def cartesian2homogeneous(cartCoords: np.ndarray) -> np.ndarray:
        homogenousCoords = np.vstack((cartCoords, np.ones(cartCoords.shape[1])))
        return homogenousCoords
    
    @staticmethod
    def homogeneous2cartesian(homogeneousCoords: np.ndarray) -> np.ndarray:
        cartCoords = homogeneousCoords[0:2, :]
        return cartCoords / homogeneousCoords[2, :]
    
    
    