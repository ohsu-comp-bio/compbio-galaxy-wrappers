'''
Keep track of how long functions take to complete. This is meant to track functions that are 
called multiple times. It keeps track of how long each call takes. And supports tracking of
multiple functions simultaneously. 

See test_benchmarking.py for example usage. 

Created on Sep. 25, 2024

@author: pleyte
'''
import logging
import statistics
from time import perf_counter

class Benchmark():
    def __init__(self, name):
        self.name = name
        self.start = None
        self.stop = None        
        self.durations = []

class Benchmarking(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)
        self._benchmarks = dict()
        self._last_name_updated = None
    
    def start(self, name):
        """
        Start timing an event
        """
        assert name, "name must be specified"
        benchmark = self._benchmarks.get(name)

        if not benchmark:
            benchmark = Benchmark(name)
            self._benchmarks[name] = benchmark 

        if benchmark.start or benchmark.stop:
            raise ValueError(f"The stopwatch for {name} needs to be stopped or cancelled before it can be started again.")

        self._last_name_updated = name

        benchmark.start = perf_counter()

    def stop(self, name):
        """
        Stop the timer for an event
        """
        benchmark = self._benchmarks.get(name)
        
        if not benchmark or not benchmark.start:
            raise ValueError(f"Unable to stop because the stopwatch for {name} not been started")
        
        benchmark.stop = perf_counter()
        
        duration = benchmark.stop - benchmark.start
        
        benchmark.durations.append(duration)
        
        benchmark.start = None
        benchmark.stop = None
        
        self._last_name_updated = None
        
    def cancel(self, name):
        """
        End an event timer without recording the duration.
        """
        benchmark = self._benchmarks.get(name)
        
        if not benchmark or not benchmark.start:
            raise ValueError(f"Unable to cancel because the stopwatch for {name} not been started") 
       
        benchmark.start = None
        benchmark.stop = None
    
    def cancel_last(self):
        """
        Cancel the stopwatch for whatever event was most recently started. 
        """
        self.cancel(self._last_name_updated)

    def clear(self):
        """
        Remove all events being timed 
        """
        self._benchmarks.clear()
    
    def get_time_total(self, name):
        """
        Return the sum of all event times in milliseconds 
        """
        benchmark = self._benchmarks.get(name) 
        return sum(benchmark.durations) * 1000.0
    
    def get_time_average(self, name):
        """
        Return the average event time in milliseconds
        """
        benchmark = self._benchmarks.get(name)
        if len(benchmark.durations) == 0:
            return 0
        
        return statistics.mean(benchmark.durations) * 1000.0
    
    def get_names(self):
        """
        Return the names of the events being monitored
        """        
        return self._benchmarks.keys()