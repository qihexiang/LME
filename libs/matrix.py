from scipy.spatial.transform import Rotation as R
import numpy as np

def mirror_matrix(law_vector):
    [x, y, z] = law_vector / np.linalg.norm(law_vector)
    return np.array(
        [
            [1 - 2 * (x**2), -2 * x * y, -2 * x * z],
            [-2 * x * y, 1 - 2 * (y**2), -2 * y * z],
            [-2 * x * z, -2 * y * z, 1 - 2 * (z**2)],
        ],
        dtype="float64",
    )


def rotate_matrix(axis, angle):
    axis = np.array(axis, dtype="float64")
    axis = axis / np.linalg.norm(axis)
    return R.from_rotvec(axis * angle, degrees=True).as_matrix()
