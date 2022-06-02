import unittest
import sys

from random import Random
random = Random(1337)

sys.path.append('.')
sys.path.append('./src')

# local files
from src.utility import *

class TestSetup(unittest.TestCase):

    def test_de2bi(self):
        return True;

    def test_get_hamiltonion(self):
        return True;

    def test_get_prob_from_hamiltonion(self):
        return True;

    def test_nchoosek(self):
        # given
        n = 1
        k = 1
        # when
        res = nchoosek(n,k)
        # then
        self.assertEqual(res,1)

        # given
        n = 10
        k = 6
        # when
        res = nchoosek(n,k)
        # then
        self.assertEqual(res,210)


    def test_gausswin(self):
        # given
        # L = random.randint(1, 20)
        L = 32
        # when
        res = gausswin(L)
        # then the result should be an array with the same size as the input length
        self.assertEqual(len(res), L)

    def test_calculate_entropy_from_prob_dist(self):
        prob_dist = np.zeros(4);
        prob_dist[0] = 1
        entropy = calculate_entropy_from_prob_dist(prob_dist)
        self.assertEqual(entropy, 0)
        
        # coin flip
        prob_dist = np.ones(2);
        prob_dist = prob_dist / np.sum(prob_dist)
        entropy = calculate_entropy_from_prob_dist(prob_dist)
        self.assertEqual(entropy, 1)
        
        # weighted coin flip
        prob_dist = np.array([0.7, 0.3])
        entropy = calculate_entropy_from_prob_dist(prob_dist)
        self.assertLess(entropy, 1)
        
        # dice roll
        prob_dist = np.ones(6);
        prob_dist = prob_dist / np.sum(prob_dist)
        entropy = calculate_entropy_from_prob_dist(prob_dist)
        self.assertAlmostEqual(entropy, 2.585, 3)

if __name__ == '__main__':
    unittest.main()