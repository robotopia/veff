test: 0437.par
	../veff -p 0437.par -e 56559.878 -f 192 -g -v -s /usr/local/share/jpl/de430.bsp

parabfit: ss.out
	../parabfit --orig=257,257 \
		        --res=0.195,0.195 \
				--xrange=-20,20 \
				--yrange=0,20 \
				--mask=1.2,1.6 \
				--omask=4.0,4.0 \
				--pdist=0.5,0.5 \
				--dB \
				--logspace \
				--curves=0.1:2.0:1000 \
				--q1 --q2 \
				--cbrange=-55,-38 \
				--out=mwa_cropped_Q12.hg.txt \
				--hggpi=mwa_cropped_Q12.hg.gpi \
				--crop=mwa_cropped.ss \
				--ssgpi=mwa_cropped.ss.gpi \
				$<

ss_test: ss_test.c ../ss.o
	gcc -o $@ $^ -lm

hg_test: hg_test.c ../ss.o ../hough.o
	gcc -o $@ $^ -lm
