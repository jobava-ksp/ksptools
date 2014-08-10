import datetime
import re

def toseconds(*args):
    return KerbalTime.from_ydhms(*args).total_seconds

class KerbalTime(object):
    
    seconds_per_year = 9203545.0
    seconds_per_day = 21600.0
    seconds_per_hour = 3600.0
    seconds_per_minute = 60.0
    seconds_per_microsecond = 1e-6  
    
    days_per_year = 426
    hours_per_day = 6
    minutes_per_hour = 60
    seconds_per_minute = 60
    
    def __init__(self, seconds):
        self.total_seconds = seconds
    
    @property
    def years(self):
        return int(self.total_seconds // KerbalTime.seconds_per_year)
        
    @property
    def days(self):
        return int(self.total_seconds // KerbalTime.seconds_per_day) % KerbalTime.days_per_year
    
    @property
    def hours(self):
        return int(self.total_seconds // KerbalTime.seconds_per_hour) % KerbalTime.hours_per_day
    
    @property
    def minutes(self):
        return int(self.total_seconds // KerbalTime.seconds_per_minute) % KerbalTime.minutes_per_hour
    
    @property
    def seconds(self):
        return int(self.total_seconds) % KerbalTime.seconds_per_minute
    
    @property
    def microsecnds(self):
        return int(self.total_seconds // seconds_per_microsecond)
    
    @classmethod
    def from_ydhms(cls, years, days=0., hours=0., minutes=0., seconds=0.):
        return cls(sum([
                years * cls.seconds_per_year,
                days * cls.seconds_per_day,
                hours * cls.seconds_per_hour,
                minutes * cls.seconds_per_minute,
                seconds]))
    
    def __add__(self, other):
        return KerbalTime(self.total_seconds + other.total_seconds)
    
    def __sub__(self, other):
        return KerbalTime(self.total_seconds - other.total_seconds)
    
    
class KerbalDate(KerbalTime):
    def __init__(self, timeexpr):
        if isinstance(timeexpr, str):
            raise NotImplementedError
        elif isinstance(timeexpr, KerbalTime):
            KerbalTime.__init__(self, timeexpr.total_seconds)
        else:
            KerbalTime.__init__(self, timeexpr)
