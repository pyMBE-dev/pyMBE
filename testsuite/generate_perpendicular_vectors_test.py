import numpy as np
import pyMBE
from itertools import combinations
pmb = pyMBE.pymbe_library()

def check_if_different_perpendicular_vectors_are_generated(vector,magnitude,n=50):
    """
    Checks if pmb.generate_trial_perpendicular_vector generates perpendicular vectors to `vector`.

    Args:
        vector(`array`,`float`): vector to which new perpendicular vectors will be generated.
        magnitude(`float`): magnitude of the perpendicular vectors.
        n(`int`): number of perpendicular vectors to be generated.
    """
    perpendicular_vectors = []
    for _ in range(n):
        perpendicular_vector = pmb.generate_trial_perpendicular_vector(vector=vector,
                                                                        magnitude=magnitude)
        perpendicular_vectors.append(perpendicular_vector)
    # Check that the function creates {n} perpendicular vectors to an input vector
    for pvector in perpendicular_vectors:
        np.testing.assert_almost_equal(actual = np.dot(pvector, vector),
                                desired = 0,
                                decimal = 5,
                                verbose = True)
    # Check that the {n} perpendicular vectors are different
    for vector_pair in combinations(perpendicular_vectors, 2):
        if np.array_equal(vector_pair[0],vector_pair[1]):
            raise Exception(f"Error: pmb.generate_trial_perpendicular_vector two equal perpendicular vectors v1 = {vector_pair[0]} v2 = {vector_pair[1]}")
    # Check that the perpendicular vectors have the same magnitude as the input magnitude
    for pvector in perpendicular_vectors:
        np.testing.assert_almost_equal(actual = np.linalg.norm(pvector), 
                                        desired = magnitude, 
                                        decimal = 5, 
                                        verbose = True)
    return

print(f"*** generate_trial_perpendicular_vector unit tests ***")
print(f"*** Unit test: Check that the function creates perpendicular vectors to an arbitrary vector of the same magnitude  ***")
vector = pmb.generate_random_points_in_a_sphere(center=[0,0,0],
                                                radius=1, 
                                                n_samples=1, 
                                                seed=None, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                    magnitude=1)
print(f"*** Unit test passed ***")
print(f"*** Unit test: Check that the function creates perpendicular vectors to a vector with an arbitrary origin ***")
vector = pmb.generate_random_points_in_a_sphere(center=[1,2,3],
                                                radius=1, 
                                                n_samples=1, 
                                                seed=None, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                    magnitude=1)
print(f"*** Unit test passed ***")
print(f"*** Unit test: Check that the function creates perpendicular vectors with a different magnitude  ***")
vector = pmb.generate_random_points_in_a_sphere(center=[1,2,3],
                                                radius=2, 
                                                n_samples=1, 
                                                seed=None, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                        magnitude=2)
print(f"*** Unit test passed ***")
print(f"*** All unit tests passed ***")