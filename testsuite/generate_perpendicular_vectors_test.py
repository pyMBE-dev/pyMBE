import numpy as np
import pyMBE
pmb = pyMBE.pymbe_library()

#Defining variables for test
origin = [1, 2, 3]
magnitude = 1
n = 10
tolerance = 1e-3
original_vector = pmb.generate_trialvectors(center=origin,
                                            radius=magnitude, 
                                            n_samples=1, 
                                            seed=None, 
                                            on_surface=True)[0]
perpendicular_vectors = []
for _ in range(n):
    perpendicular_vector = pmb.generate_trial_perpendicular_vector(vector=original_vector,
                                                                    origin=origin, 
                                                                    magnitude=magnitude)
    perpendicular_vectors.append(perpendicular_vector)

print(f"*** generate_trial_perpendicular_vector unit tests ***")
print(f"*** Unit test: Check that the function creates {n} perpendicular vectors to an input vector. ***")
for vector in perpendicular_vectors:
    np.testing.assert_almost_equal(actual = np.dot(original_vector, vector),
                            desired = 0,
                            decimal = 5,
                            verbose = True)

print(f"***Unit test passed***")
print(f"***Unit test: Check that the {n} perpendicular vectors are different.***")
unique_set = set()
for vector in perpendicular_vectors:
    vector_tuple = tuple(round(coord, 2) for coord in vector)
    unique_set.add(vector_tuple)
np.testing.assert_equal(actual = len(unique_set), 
                        desired = n, 
                        verbose = True)
print(f"***Unit test passed***")
   
print(f"***Unit test: Check that the perpendicular vectors have the same magnitude as the input magnitude.***")
for vector in perpendicular_vectors:
    np.testing.assert_almost_equal(actual = np.linalg.norm(vector), 
                                    desired = magnitude, 
                                    decimal = 5, 
                                    verbose = True)
print(f"***Unit test passed")
print(f"*** All unit tests passed ***")