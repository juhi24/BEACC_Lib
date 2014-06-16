import time
import h5py
import read	
	

def great_hdf5(date):
	date_str = time.strftime("%Y%m%d",date)
	f = h5py.File(date_str + ".hdf5", "w")
	dset = f.create_dataset("basetime", (1,), dtype='i')
	#read date from different device
	dev = [read.hotplate,read.pluvio,read.jenoptik]
	print dev
	for i in dev:
		data = i(date)
		names = list(data.columns.values)
		for name in names:
			data_values = data[name].values
			data_values = data_values.tolist()
			#data_values = [float(i) for i in data_values]
			#print type(data_values)
			dset = f.create_dataset(name, data=data_values)

	
def main():
	tmp=(2014,02,07,0,0,0,0,0,0)
	date=time.mktime(tmp)
	date = time.gmtime(date)
	great_hdf5(date)

main()
