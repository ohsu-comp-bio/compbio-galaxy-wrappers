'''
Test the Benchmarking class

Created on Sep. 25, 2024

@author: pleyte
'''
import math
import time
import unittest

from edu.ohsu.compbio.txeff.util.benchmarking import Benchmarking


class TestBenchmarking(unittest.TestCase):
    '''
    Test the Benchmarking class
    '''    
    def test__example1(self):
        benchmarking = Benchmarking()
                
        n = 10
        
        # Number of seconds that each process will sleep for
        # A's timer runs while it sleeps and its timer continues to run while B sleeps.    
        a_sleep = 1/4
        b_sleep = 1/8
        
        for x in range(n):
            benchmarking.start("A")
            time.sleep(a_sleep)
        
            benchmarking.start("B")
            time.sleep(b_sleep)
            benchmarking.stop("B")
            
            benchmarking.stop("A")

        # Calculate expected amount of time spent sleeping and convert to milliseconds
        a_expected_total = ( (n * a_sleep) + (n * b_sleep) ) * 1000.0
        a_expected_average = a_expected_total / n
        a_actual_total = benchmarking.get_time_total('A')
        a_actual_average = benchmarking.get_time_average('A')
        
        # Since time spent sleeping isn't exact, compare total milliseconds with a tolerance of +/- 100.0ms     
        self.assertTrue(math.isclose(a_expected_total, a_actual_total, abs_tol=1e+2), f'Total times are not close: {a_expected_total} and {a_actual_total}')
        
        # Average is smaller so compare average milliseconds with a tolerance of +/- 10.0ms
        self.assertTrue(math.isclose(a_expected_average, a_actual_average, abs_tol=1e+1), f'Total times are not close: {a_expected_average} and {a_actual_average}')
        
        b_expected_total = n * b_sleep * 1000.0
        b_expected_average = b_expected_total / n 
        b_actual_total = benchmarking.get_time_total('B')
        b_actual_average = benchmarking.get_time_average('B')

        self.assertTrue(math.isclose(b_expected_total, b_actual_total, abs_tol=1e+2), f'Total times are not close: {b_expected_total} and {b_actual_total}')        
        self.assertTrue(math.isclose(b_expected_average, b_actual_average, abs_tol=1e+1), f'Total times are not close: {b_expected_average} and {b_actual_average}')