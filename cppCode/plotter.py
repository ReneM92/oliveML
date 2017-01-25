from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import glob
import math


def plotAll(fname, xMin, xMax, yMin, yMax, dt, simLength, delimiter = '\t', skiprows = 1):
	with open(fname, 'r') as f:
		firstLine = f.readline().strip()
	headers = firstLine.split()
	nVariables = len(headers)
	
	totalSteps = int(simLength / dt)
	xStart = int(xMin / dt)
	xEnd = int(xMax / dt)

	times = np.arange(0, simLength, dt)

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	values = np.loadtxt(fname = fname, delimiter = delimiter, skiprows = skiprows);
	#print(values)
	for i in range(1,nVariables):
		#print(values[xStart: xEnd, i])
		ax1.plot(times[xStart: xEnd], values[xStart: xEnd, i], label = headers[i])

	ax1.set_xlabel('time(s)')
	ax1.set_ylabel('value')
	ax1.legend(loc = 'best')
	ax1.grid(True)
	ax1.set_xlim([xMin, xMax])
	ax1.set_ylim([yMin, yMax])

	plt.savefig('./images/results.png')
	plt.close()

def plotAllTauAndInf(delimiter = '\t', skiprows = 1):
	"""
		dendrite
	"""
	dendTau_k = np.loadtxt(fname = 'channel_summary/dend_tau_k.dat', delimiter = delimiter, skiprows = skiprows)
	dendInf_k = np.loadtxt(fname = 'channel_summary/dend_inf_k.dat', delimiter = delimiter, skiprows = skiprows)
	dendTau_l = np.loadtxt(fname = 'channel_summary/dend_tau_l.dat', delimiter = delimiter, skiprows = skiprows)
	dendInf_l = np.loadtxt(fname = 'channel_summary/dend_inf_l.dat', delimiter = delimiter, skiprows = skiprows)

	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111)
	ax1.plot(dendInf_k[:,0], dendInf_k[:,1], label = "k")
	ax1.plot(dendInf_l[:,0], dendInf_l[:,1], label = 'l')
	ax1.set_xlabel("voltage")
	ax1.set_ylabel("inf")
	ax1.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax1.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_inf_cal.png')

	fig2 = plt.figure(2)
	ax2 = fig2.add_subplot(111)
	ax2.plot(dendTau_k[:,0], dendTau_k[:,1], label = "k")
	ax2.plot(dendTau_l[:,0], dendTau_l[:,1], label = 'l')
	ax2.set_xlabel("voltage")
	ax2.set_ylabel("tau")
	ax2.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax2.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_tau_cal.png')

	dendTau_r = np.loadtxt(fname = 'channel_summary/dend_tau_r.dat', delimiter = delimiter, skiprows = skiprows)
	dendInf_r = np.loadtxt(fname = 'channel_summary/dend_inf_r.dat', delimiter = delimiter, skiprows = skiprows)

	fig3 = plt.figure(3)
	ax3 = fig3.add_subplot(111)
	ax3.plot(dendInf_r[:,0], dendInf_r[:,1], label = 'r')
	ax3.set_xlabel("voltage")
	ax3.set_ylabel("inf")
	ax3.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax3.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_inf_cah.png')

	fig4 = plt.figure(4)
	ax4 = fig4.add_subplot(111)
	ax4.plot(dendTau_r[:,0], dendTau_r[:,1], label = 'r')
	ax4.set_xlabel("voltage")
	ax4.set_ylabel("tau")
	ax4.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax4.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_tau_cah.png')

	dendTau_n = np.loadtxt(fname = 'channel_summary/dend_tau_n.dat', delimiter = delimiter, skiprows = skiprows)
	dendInf_n = np.loadtxt(fname = 'channel_summary/dend_inf_n.dat', delimiter = delimiter, skiprows = skiprows)

	fig5 = plt.figure(5)
	ax5 = fig5.add_subplot(111)
	ax5.plot(dendInf_n[:,0], dendInf_n[:,1], label = 'n')
	ax5.set_xlabel("voltage")
	ax5.set_ylabel("inf")
	ax5.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax5.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_inf_h.png')

	fig6 = plt.figure(6)
	ax6 = fig6.add_subplot(111)
	ax6.plot(dendTau_n[:,0], dendTau_n[:,1], label = 'n')
	ax6.set_xlabel("voltage")
	ax6.set_ylabel("tau")
	ax6.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax6.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./dend_tau_h.png')

	"""
		axon
	"""
	axonInf_m = np.loadtxt(fname = 'channel_summary/axon_inf_m.dat', delimiter = delimiter, skiprows = skiprows)
	axonTau_h = np.loadtxt(fname = 'channel_summary/axon_tau_h.dat', delimiter = delimiter, skiprows = skiprows)
	axonInf_h = np.loadtxt(fname = 'channel_summary/axon_inf_h.dat', delimiter = delimiter, skiprows = skiprows)

	fig7 = plt.figure(7)
	ax7 = fig7.add_subplot(111)
	ax7.plot(axonInf_m[:,0], axonInf_m[:,1], label = "m")
	ax7.plot(axonInf_h[:,0], axonInf_h[:,1], label = 'h')
	ax7.set_xlabel("voltage")
	ax7.set_ylabel("inf")
	ax7.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax7.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./axon_inf_na_a.png')

	fig8 = plt.figure(8)
	ax8 = fig8.add_subplot(111)
	ax8.plot(axonTau_h[:,0], axonTau_h[:,1], label = 'h')
	ax8.set_xlabel("voltage")
	ax8.set_ylabel("tau")
	ax8.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax8.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./axon_tau_na_a.png')

	axonTau_nk = np.loadtxt(fname = 'channel_summary/axon_tau_nk.dat', delimiter = delimiter, skiprows = skiprows)
	axonInf_nk = np.loadtxt(fname = 'channel_summary/axon_inf_nk.dat', delimiter = delimiter, skiprows = skiprows)

	fig9 = plt.figure(9)
	ax9 = fig9.add_subplot(111)
	ax9.plot(axonInf_nk[:,0], axonInf_nk[:,1], label = 'nk')
	ax9.set_xlabel("voltage")
	ax9.set_ylabel("inf")
	ax9.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax9.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./axon_inf_k.png')

	fig10 = plt.figure(10)
	ax10 = fig10.add_subplot(111)
	ax10.plot(axonTau_nk[:,0], axonTau_nk[:,1], label = 'nk')
	ax10.set_xlabel("voltage")
	ax10.set_ylabel("tau")
	ax10.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax10.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./axon_tau_k.png')

	"""
		soma
	"""
	somaInf_m = np.loadtxt(fname = 'channel_summary/soma_inf_m.dat', delimiter = delimiter, skiprows = skiprows)
	somaTau_h = np.loadtxt(fname = 'channel_summary/soma_tau_h.dat', delimiter = delimiter, skiprows = skiprows)
	somaInf_h = np.loadtxt(fname = 'channel_summary/soma_inf_h.dat', delimiter = delimiter, skiprows = skiprows)

	fig11 = plt.figure(11)
	ax11 = fig11.add_subplot(111)
	ax11.plot(somaInf_m[:,0], somaInf_m[:,1], label = "m")
	ax11.plot(somaInf_h[:,0], somaInf_h[:,1], label = 'h')
	ax11.set_xlabel("voltage")
	ax11.set_ylabel("inf")
	ax11.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax11.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_inf_na_s.png')

	fig12 = plt.figure(12)
	ax12 = fig12.add_subplot(111)
	ax12.plot(somaTau_h[:,0], somaTau_h[:,1], label = 'h')
	ax12.set_xlabel("voltage")
	ax12.set_ylabel("tau")
	ax12.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax12.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_tau_na_s.png')


	somaTau_n = np.loadtxt(fname = 'channel_summary/soma_tau_n.dat', delimiter = delimiter, skiprows = skiprows)
	somaInf_n = np.loadtxt(fname = 'channel_summary/soma_inf_n.dat', delimiter = delimiter, skiprows = skiprows)

	fig13 = plt.figure(13)
	ax13 = fig13.add_subplot(111)
	ax13.plot(somaInf_n[:,0], somaInf_n[:,1], label = "n")
	ax13.set_xlabel("voltage")
	ax13.set_ylabel("inf")
	ax13.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax13.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_inf_kdr.png')

	fig14 = plt.figure(14)
	ax14 = fig14.add_subplot(111)
	ax14.plot(somaTau_h[:,0], somaTau_h[:,1], label = 'n')
	ax14.set_xlabel("voltage")
	ax14.set_ylabel("tau")
	ax14.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax14.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_tau_kdr.png')

	somaTau_nk = np.loadtxt(fname = 'channel_summary/soma_tau_nk.dat', delimiter = delimiter, skiprows = skiprows)
	somaInf_nk = np.loadtxt(fname = 'channel_summary/soma_inf_nk.dat', delimiter = delimiter, skiprows = skiprows)

	fig15 = plt.figure(15)
	ax15 = fig15.add_subplot(111)
	ax15.plot(somaInf_n[:,0], somaInf_n[:,1], label = "nk")
	ax15.set_xlabel("voltage")
	ax15.set_ylabel("inf")
	ax15.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax15.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_inf_k.png')

	fig16 = plt.figure(16)
	ax16 = fig16.add_subplot(111)
	ax16.plot(somaTau_h[:,0], somaTau_h[:,1], label = 'nk')
	ax16.set_xlabel("voltage")
	ax16.set_ylabel("tau")
	ax16.set_xlim([-110 * math.pow(10, -3), 110 * math.pow(10, -3)])
	ax16.grid(True)
	#print(dendInf_k[:,1])
	plt.savefig('./soma_tau_k.png')


def main():
	
	fname = "test.dat"
	xMin = 0
	xMax = 1000
	yMin = -100 * math.pow(10,-3)
	yMax = 100 * math.pow(10,-3)

	dt = 0.025
	simLength = 2000

	#for time in np.arange(0,0.75, 0.25):
	#	print(time)
	plotAll(fname, xMin, xMax, yMin, yMax, dt, simLength)
	#plotAllTauAndInf()


if __name__ == "__main__":
  main()
