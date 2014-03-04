import sys, os, time
from scipy import ndimage
from scipy import misc
import scipy.io
from matplotlib.image import AxesImage
from pylab import *
import numpy as np
import subprocess
import glob
import shutil
from PIL import Image, ImageDraw, ImageChops, ImageFilter
import matplotlib.pyplot as plt
import datetime,time
from os.path import join, getsize
import re
import h5py
import pandas as pd

def hotplate(date):
	date_str = time.strftime("%Y%m%d",date)
	files = glob.glob("data/Hotplate/hot_plate_100901_"+date_str+"*")
	acc = []
	time_ = []
	for filename in files:
		file_ = open(filename)
		lines = file_.readlines()
		file_.close
		for i in lines:
			var = i.split(',')
			if len(var) > 0:
				time_tmp = time.strptime(var[0],'%Y%m%d%H%M%S')
				time_tmp = time.mktime(time_tmp)
				time_.append(time_tmp)
				acc.append(float(var[-1]))
	
	#print acc
	d = {'hotplate_time' : time_, 'hotplate_accumulation': acc}
	return pd.DataFrame(d)

def jeoptic(date):
	date_str = time.strftime("%Y%m%d",date)
	files = glob.glob("data/Jenoptik/"+date_str+"*")
	snow = []
	time_ = []
	data_=[]
	for filename in files:
		file_ = open(filename)
		lines = file_.readlines()
		file_.close
		for i in lines:
			var = i.split(',')
			if len(var) > 0:
				time_tmp = time.strptime(var[0],'%Y-%m-%d %H:%M:%S')
				time_tmp = time.mktime(time_tmp)
				time_.append(time_tmp)
				snow.append(float(var[2])-0.034)
	#print acc
	d = {'jenoptik_time' : time_, 'jenoptik_snow_depth': snow}
	return pd.DataFrame(d)

def parsivel23(date):
	date_str = time.strftime("%Y%m%d",date)
	files = glob.glob("data/Parsivel23/"+date_str+"*")	
	
def pip(filename):
	d = pd.read_csv(filename, delim_whitespace=True, skiprows=8, header=3)
	return d
