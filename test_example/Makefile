test: 0437.par
	../veff -p 0437.par -e 56559.878 -f 192 -g -v -s /usr/local/share/jpl/de430.bsp

parabfit: ss.out
	../parabfit --ssdata=$< --orig=257,257 --res=0.195,0.195 --xrange=-20,20 --yrange=0,20 --mask=20,7,8 --tmax=5 --dB

