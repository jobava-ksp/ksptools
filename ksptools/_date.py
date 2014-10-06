import re
import numpy as np

time_to_sec_conv = np.array([426*21600, 21600, 3600, 60, 1])
sec_to_time_conv = np.array([1/(426*21600.0), 1/21600.0, 1/3600.0, 1/60.0, 1])

def timetosec(ktime_string):
    """
    Convert a Kerbal time string into seconds.
    :type ktime_string: str
    :rtype: float
    """
    date_regex = r'\+Y(?P<y>[0-9]+)\s*D(?P<d>[0-9]+)\s*,?\s*(?P<h>[0-9]+):(?P<m>[0-9]+):(?P<s>[0-9]+)'
    m = re.match(date_regex, ktime_string)
    time_array = np.array([float(m.group('y')), float(m.group('d')), float(m.group('h')), float(m.group('m')), float(m.group('s'))])
    return np.dot(time_to_sec_conv, time_array)


def datetosec(kdate_string):
    """
    Convert a Kerbal date string into seconds.
    :type kdate_string: str
    :rtype: float
    """
    date_regex = r'\Y(?P<y>[0-9]+)\s*D(?P<d>[0-9]+)\s*,?\s*(?P<h>[0-9]+):(?P<m>[0-9]+):(?P<s>[0-9]+)'
    m = re.match(date_regex, kdate_string)
    time_array = np.array([float(m.group('y'))-1, float(m.group('d'))-1, float(m.group('h')), float(m.group('m')), float(m.group('s'))])
    return np.dot(time_to_sec_conv, time_array)


def sectotime(time):
    """
    Convert seconds into a Kerbal time string.
    :type time: float
    :rtype: str
    """
    time_arr = np.array([int(time*sec_to_time_conv[0]), 0, 0, 0, 0])
    for i in range(1,5):
        time_arr[i] = int((time - np.dot(time_to_sec_conv[0:i], time_arr[0:i]))*sec_to_time_conv[i])
    return '+Y{} D{}, {:0>2d}:{:0>2d}:{:0>2d}'.format(*time_arr)


def sectodate(time):
    """
    Convert seconds into a Kerbal date string.
    :type time: float
    :rtype: str
    """
    time_arr = np.array([int(time*sec_to_time_conv[0]), 0, 0, 0, 0])
    for i in range(1,5):
        time_arr[i] = int((time - np.dot(time_to_sec_conv[0:i], time_arr[0:i]))*sec_to_time_conv[i])
    y,d,h,m,s = time_arr
    return 'Y{} D{}, {:0>2d}:{:0>2d}:{:0>2d}'.format(y+1,d+1,h,m,s)


