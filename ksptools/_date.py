import re
import numpy as np

timeconv = np.array([9203545, 21600, 3600, 60, 1])

def date_to_time(kdate_string):
    date_regex = r'Y(?P<y>[0-9]+)\s*D(?P<d>[0-9]+)\s*,?\s*(?P<h>[0-9]+):(?P<m>[0-9]+):(?P<s>[0-9]+)'
    m = re.match(date_regex, kdate_string)
    time_array = np.array([float(m.group('y'))-1, float(m.group('d'))-1, float(m.group('h')), float(m.group('m')), float(m.group('s'))])
    return np.dot(timeconv, time_array)


def time_to_date(time):
    y = int(time/timeconv[0])
    d = int((time - np.dot(timeconv[0:1], np.array([y])))/timeconv[1])
    h = int((time - np.dot(timeconv[0:2], np.array([y,d])))/timeconv[2])
    m = int((time - np.dot(timeconv[0:3], np.array([y,d,h])))/timeconv[3])
    s = int((time - np.dot(timeconv[0:4], np.array([y,d,h,m])))/timeconv[4])
    return 'Y{} D{}, {}:{}:{}'.format(y+1,d+1,h,m,s)


